function [art_suspects, FD_run, delta_run] = art_global_qa(Images, RealignmentFile, HeadMaskType, RepairType, Percent_thresh, mv_thresh)
% FORMAT art_global                         (v.2.6)
%
%     Art_global allows visual inspection of average intensity and
% scan to scan motion of fMRI data, and offers methods to repair outliers
% in the data. Outliers are scans whose global intensity is very different
% from the mean, or whose scan-to-scan motion is large. Thresholds that
% define the outliers are shown, and they can be adjusted by the user.
% Margins are defined around the outliers to assure that the repaired
% data will satisfy the "slowly-varying" background assumption of GLM
% models. Outliers can also be defined by user manual edits. When all
% repair parameters are set, the user writes new repaired files by using
% the Repair button in the GUI. Repairs can be done by interpolation
% between the nearest non-repaired scans (RECOMMNEDED), 
% or despike interpolation using the
% immediate before and after scans, or inserting the mean image in place
% of the image to be repaired. Repairs will change the scans marked by
% red vertical lines. Scans marked by green vertical lines are unchanged,
% but will be added to the deweighting text file.
%     Repaired images have the prefix "v" and are written to the same
% folder as the input images. The input images are preserved in place.
% Unchanged images are copied to a "v" named file, so SPM can run with
% the "v" images. The program writes a file art_repaired.txt with a list
% of all repaired images, and art_deweighted.txt with a list of all images
% to be deweighted during estimation. Deweighting is recommended for both 
% repaired images and margin images. Repairs bias contrasts lower, but 
% informally, the effect seems small when fewer than 5-10% of scans are 
% repaired. For deeper repairs, deweighting must be implemented in a
% batch script, e.g. as in the art_redo function.
%     A set of default thresholds is suggested ( 1.3% variation in global
% intensity, 0.5 mm/TR scan-to-scan motion.) Informally, these values
% are OK for small to moderate repairs. The thresholds can be reduced for
% good data (motion -> 0.3) , and should be raised for severely noisy data
% (motion -> 1.0). The values of the intensity and motion thresholds are
% linked and will change together. As a default, all images with total
% movement > 3 mm will also be marked for repair. For special situations,
% the outlier edit feature can mark additional scans. When motion 
% regressors will be used, suggest setting the motion threshold to 0.5
% and not applying the margins.
%
% For MULTIPLE SESSIONS, we suggest realigning each session separately,
% and repairing each session separately. This approach agrees more
% with SPM standard practice, and clinical subjects often move between
% sessions. This approach differs from previous versions of this program. 
%
% For batch scripts, use
% FORMAT art_global(Images, RealignmentFile, HeadMaskType, RepairType)
%    Images  = Full path name of images to be repaired.
%       For multiple sessions, all Images are in one array, e.g.SPM.xY.P
%    RealignmentFile = Full path name of realignment file
%       For multiple sessions, cell array with one name per session.
%    HeadMaskType  = 1 for SPM mask, = 4 for Automask
%    RepairType = 1 for ArtifactRepair alone (0.5 movement and add margin).
%               = 2 for Movement Adjusted images  (0.5 movement, no margin)
%               = 0 No repairs are done, bad scans are found.
%                   Listed in art_suspects.txt for motion adjustment.
%    Hardcoded actions:
%       Does repair, does not force repair of first scan.
% ----------------------------------------------------------------------

% MemoLab QA toolbox changes, August 2016
%   - additional inputs for Percent_thresh and mv_thresh
%   - specified outputs
%   - For all MemoLab changes, do a Find for "memolab"

% ----------------------------------------------------------------------
% v2.5, May 2009 pkm  - adds SPM8, RepairType=0.
% v2.6, Dec 2014 supports SPM12.

% v2.4, Mar 2009 pkm
%    Only one session allowed. Realign and repair each session separately.
%    Increment global intensity at 0.05%, instead of 0.1 in z-score.
%    Movement adaptive threshold now affects global threshold.
% v2.3, July 2008 PKM
%    Allows automatic adaptive threshold for movement
%    Compatible SPM5 and SPM2
%    Lowered default thresholds, to catch more deep breath artifacts.
% v2.2, July 2007 PKM
%    Also marks scan before a big movement for mv_out repair.
%    Allows two kinds of RepairType in batch, with and without margin.
%    Fixed logic for user masked mean
%    Adds user GUI option to repair all scans with movement > 3 mm.
%    small fix to subplots
% v2.1, Jan. 2007. PKM.
%    Prints copy of art_global figure in batch mode.
%    Allows multiple sessions in batch mode.
%       Images is full path name of all sessions, e.g. SPM.xY.P
%       Realignment file is one cell per session with full rp file name.
%    For SPM5, finds size of VY.dim.  ( Compatible with SPM2)
%
% Paul Mazaika, September 2006.

% -----------------------
% Initialize, begin loop
% -----------------------

pfig = [];
% Configure while preserving old SPM versions
spmv = spm('Ver'); spm_ver = 'spm5';  % chooses spm_select to read vols
if (strcmp(spmv,'SPM2')) spm_ver = 'spm2'; end
if (strcmp(spmv,'SPM2') || strcmp(spmv,'SPM5')) spm_defaults;
    else spm('Defaults','fmri'); end

% % ------------------------
% % Default values for outliers     % commented by MemoLab to take input thresholds instead
% % ------------------------
% % When std is very small, set a minimum threshold based on expected physiological
% % noise. Scanners have about 0.1% error when fully in spec. 
% % Gray matter physiology has ~ 1% range, ~0.5% RMS variation from mean. 
% % For 500 samples, expect a 3-sigma case, so values over 1.5% are
% % suspicious as non-physiological noise. Data within that range are not
% % outliers. Set the default minimum percent variation to be suspicious...
%       Percent_thresh = 1.3; 
% %  Alternatively, deviations over 2*std are outliers, if std is not very small.
%       z_thresh = 2;  % Currently not used for default.
% % Large intravolume motion may cause image reconstruction
% % errors, and fast motion may cause spin history effects.
% % Guess at allowable motion within a TR. For good subjects,
% % would like this value to be as low as 0.3. For clinical subjects,
% % set this threshold higher.
%       mv_thresh = 0.5;  % try 0.3 for subjects with intervals with low noise
%                         % try 1.0 for severely noisy subjects   

% ------------------------
% Collect files
% ------------------------
if nargin > 0
    % num_sess = 1;   % Only one session;
    num_sess = size(RealignmentFile,1);
    global_type_flag = HeadMaskType;
    realignfile = 1;
    P{1} = Images;
    M = [];
    for i = 1:num_sess
        %P{i} = Images{i};
        mvmt_file = char(RealignmentFile(i,:));
        M{i} = load(mvmt_file);
        %[mv_path,mv_name,mv_ext] = fileparts(mvmt_file);
        % M{1} = mv_data;
    end
    %repair1_flag = 0;   % Only repair scan 1 when necessary
    repair1_flag = 1;   % Force repair scan 1  (Reisslab)
    GoRepair = 1;       % Automatic Repair
    if nargin >= 4      % To stay backward compatible with v2.1 - set to >= instead of == for extra inputs (MemoLab)
        if RepairType == 1
            GoRepair = 1;
        elseif RepairType == 2
%             mv_thresh = 0.5; % never overwrite input thresholds (MemoLab)
            GoRepair = 2;
        elseif RepairType == 0  % no repairs done, bad scans found.
%             mv_thresh = 0.5; % never overwrite input thresholds (MemoLab)
            GoRepair = 4;
        end
    end
else
    %num_sess = spm_input('How many sessions?',1,'n',1,1);
    num_sess = 1;
    global_type_flag = spm_input('Which global mean to use?', 1, 'm', ...
    'Regular SPM mask  | Auto ( Generates ArtifactMask and can Calculate Movement )| User Mask  ( in same voxel space as functional images ) | Every Voxel',...
              [1 4 3 2], 2);
    if global_type_flag==3
         if strcmp(spm_ver,'spm5')
            mask = spm_select(1, 'image', 'Select mask image in functional space'); 
         else  % spm2 version
            mask = spm_get(1, '.img', 'Select mask image in functional space');
         end
    end
    % If there are no realignment files available, compute some instead.
    realignfile = 1;   % Default case is like artdetect4.
    if global_type_flag == 4
        realignfile = spm_input('Have realignment files?',1, 'b', ' Yes | No ', [ 1 0 ], 1);
    end

    M = []; P = [];
    if strcmp(spm_ver,'spm5')
        for i = 1:num_sess
            P{i} = spm_select(Inf,'image',['Select data images for session'  num2str(i) ':']);
            if realignfile == 1
                mvmt_file = spm_select(1,'any',['Select movement params file for session' num2str(i) ':']);      
                M{i} = load(mvmt_file);
            end
        end
    else   % spm2 version
        for i = 1:num_sess
            P{i} = spm_get(Inf,'.img',['Select data images for session'  num2str(i) ':']);
            if realignfile == 1
                mvmt_file = spm_get(1,'.txt',['Select movement params file for session' num2str(i) ':']);
                M{i} = load(mvmt_file);
            end
        end
    end
    repair1_flag = spm_input('Always repair 1st scan of each session?', '+1', 'y/n', [1 0], 1);
    GoRepair = 0;       % GUI repair
    RepairType = 2;     % Allows adaptive threshold
end  % End of GUI loop

if global_type_flag==3
    maskY = spm_read_vols(spm_vol(mask));
    %maskXYZmm = maskXYZmm(:,find(maskY==max(max(max(maskY)))));
    maskcount = sum(sum(sum(maskY)));  %  Number of voxels in mask.
    voxelcount = prod(size(maskY));    %  Number of voxels in 3D volume.
end
if global_type_flag == 4   %  Automask option
    disp('Generated mask image is written to file ArtifactMask.img.')
    Pnames = P{1};
    Automask = art_automask(Pnames(1,:),-1,1);
    maskcount = sum(sum(sum(Automask)));  %  Number of voxels in mask.
    voxelcount = prod(size(Automask));    %  Number of voxels in 3D volume.
end
spm_input('!DeleteInputObj');

P = char(P);
%if nargin == 0
    mv_data = [];
    for i = 1:length(M)
        mv_data = vertcat(mv_data,M{i});
    end
%end


% -------------------------
% get file identifiers and Global values
% -------------------------

fprintf('%-4s: ','Mapping files...')                                  
VY     = spm_vol(P);
fprintf('%3s\n','...done')                                          

temp = any(diff(cat(1,VY.dim),1,1),1);
if strcmp(spm_ver,'spm5')
     % or could test length(temp) == 3
     if ~isempty(find(diff(cat(1,VY.dim)) ~= 0 ))   
 	    error('images do not all have the same dimensions (SPM5)')
     end
elseif length(temp) == 4       % SPM2 case
     if any(any(diff(cat(1,VY.dim),1,1),1)&[1,1,1,0])
         error('images do not all have the same dimensions')
     end
end

nscans = size(P,1);
%keyboard;
% ------------------------
% Compute Global variate
%--------------------------

%GM     = 100;
g      = zeros(nscans,1);

fprintf('%-4s: %3s','Calculating globals...',' ')
if global_type_flag==1  % regular mean
    for i  = 1:nscans  
	    g(i) = spm_global(VY(i));
    end
elseif global_type_flag==2  % every voxel
    for i = 1:nscans
        g(i) = mean(mean(mean(spm_read_vols(VY(i)))));
    end
elseif global_type_flag == 3 % user masked mean
     Y = spm_read_vols(VY(1));
    %voxelcount = prod(size(Y));
    %vinv = inv(VY(1).mat);
    %[dummy, idx_to_mask] = intersect(XYZmm', maskXYZmm', 'rows');
    %maskcount = length(idx_to_mask);
    for i = 1:nscans
        Y = spm_read_vols(VY(i)); 
        Y = Y.*maskY;
        %Y(idx_to_mask) = [];   
        g(i) = mean(mean(mean(Y)))*voxelcount/maskcount;
        %g(i) = mean(Y(idx_to_mask));
    end
else   %  global_type_flag == 4  %  auto mask
    for i = 1:nscans
        Y = spm_read_vols(VY(i));
        Y = Y.*Automask;
        if realignfile == 0
            output = art_centroid(Y);
            centroiddata(i,1:3) = output(2:4);
            g(i) = output(1)*voxelcount/maskcount;
        else     % realignfile == 1
            g(i) = mean(mean(mean(Y)))*voxelcount/maskcount;
        end
    end
    % If computing approximate translation alignment on the fly...
    %   centroid was computed in voxels
    %   voxel size is VY(1).mat(1,1), (2,2), (3,3).
    %   calculate distance from mean as our realignment estimate
    %   set rotation parameters to zero.
    if realignfile == 0    % change to error values instead of means.
        centroidmean = mean(centroiddata,1);
        for i = 1:nscans
            mv0data(i,:) = - centroiddata(i,:) + centroidmean;
        end
        % THIS MAY FLIP L-R  (x translation)
        mv_data(1:nscans,1) = mv0data(1:nscans,1)*VY(1).mat(1,1);
        mv_data(1:nscans,2) = mv0data(1:nscans,2)*VY(1).mat(2,2);
        mv_data(1:nscans,3) = mv0data(1:nscans,3)*VY(1).mat(3,3);
        mv_data(1:nscans,4:6) = 0;
    end
end

% Convert rotation movement to degrees
mv_data(:,4:6)= mv_data(:,4:6)*180/pi; 
    
    
fprintf('%s%3s\n','...done\n')
if global_type_flag==3
    fprintf('\n%g voxels were in the user mask.\n', maskcount)
end
if global_type_flag==4
    fprintf('\n%g voxels were in the auto generated mask.\n', maskcount)
end

% ------------------------
% Compute default out indices by z-score, or by Percent-level is std is small.
% ------------------------ 
%  Consider values > Percent_thresh as outliers (instead of z_thresh*gsigma) if std is small.
    gsigma = std(g);
    gmean = mean(g);
    pctmap = 100*gsigma/gmean;
    mincount = Percent_thresh*gmean/100;
    %z_thresh = max( z_thresh, mincount/gsigma );
    z_thresh = mincount/gsigma;        % Default value is PercentThresh.
    z_thresh = 0.1*round(z_thresh*10); % Round to nearest 0.1 Z-score value
    zscoreA = (g - mean(g))./std(g);  % in case Matlab zscore is not available
    glout_idx = (find(abs(zscoreA) > z_thresh))';

% ------------------------
% Compute default out indices from rapid movement
% ------------------------ 
%   % Rotation measure assumes voxel is 65 mm from origin of rotation.
    if realignfile == 1 | realignfile == 0
        delta = zeros(nscans,1);  % Mean square displacement in two scans
        for i = 2:nscans
            delta(i,1) = (mv_data(i-1,1) - mv_data(i,1))^2 +...
                    (mv_data(i-1,2) - mv_data(i,2))^2 +...
                    (mv_data(i-1,3) - mv_data(i,3))^2 +...
                    1.28*(mv_data(i-1,4) - mv_data(i,4))^2 +...
                    1.28*(mv_data(i-1,5) - mv_data(i,5))^2 +...
                    1.28*(mv_data(i-1,6) - mv_data(i,6))^2;
            delta(i,1) = sqrt(delta(i,1));
        end
    end
    
     % Also name the scans before the big motions (v2.2 fix)
    deltaw = zeros(nscans,1);
    for i = 1:nscans-1
        deltaw(i) = max(delta(i), delta(i+1));
    end
    delta(1:nscans-1,1) = deltaw(1:nscans-1,1);
    
%     % Adapt the threshold  (v2.3 fix)         - commented out for MemoLab QA
%     if RepairType == 2 | GoRepair == 4
%         delsort = sort(delta);
%         if delsort(round(0.75*nscans)) > mv_thresh
%             mv_thresh = min(1.0,delsort(round(0.75*nscans)));
%             words = ['Automatic adjustment of movement threshold to ' num2str(mv_thresh)];
%             disp(words)
%             Percent_thresh = mv_thresh + 0.8;    % v2.4
%         end
%     end
    
    mvout_idx = find(delta > mv_thresh)';
    
    % Total repair list
    out_idx = unique([mvout_idx glout_idx]);
    if repair1_flag == 1
        out_idx = unique([ 1 out_idx]);
    end
    % Initial deweight list before margins
    outdw_idx = out_idx; 
    % Initial clip list without removing large displacements
    clipout_idx = out_idx;
     

% -----------------------
% Draw initial figure
% -----------------------


figure('Units', 'normalized', 'Position', [0.2 0.2 0.6 0.7]);
rng = max(g) - min(g);   % was range(g);
pfig = gcf;
% Don't show figure in batch runs
if (nargin > 0); set(pfig,'Visible','off'); end

subplot(5,1,1);
plot(g);
%xlabel(['artifact index list [' int2str(out_idx') ']'], 'FontSize', 8, 'Color','r');
%ylabel(['Range = ' num2str(rng)], 'FontSize', 8);
ylabel('Global Avg. Signal');
xlabel('Red vertical lines are to depair. Green vertical lines are to deweight.');
title('ArtifactRepair GUI to repair outliers and identify scans to deweight');
%if ( global_type_flag == 1 ) title('Global Mean - Regular SPM'); end
%if ( global_type_flag == 2 ) title('Global Mean - Every Voxel'); end
%if ( global_type_flag == 3 ) title('Global Mean - User Defined Mask'); end
%if ( global_type_flag == 4 ) title('Global Mean - Generated ArtifactMask'); end

% Add vertical exclusion lines to the global intensity plot
axes_lim = get(gca, 'YLim');
axes_height = [axes_lim(1) axes_lim(2)];
for i = 1:length(outdw_idx)   % Scans to be Deweighted
    line((outdw_idx(i)*ones(1, 2)), axes_height, 'Color', 'g');
end
if GoRepair == 2
    for i = 1:length(outdw_idx)   % Scans to be Deweighted
        line((outdw_idx(i)*ones(1, 2)), axes_height, 'Color', 'r');
    end
end
subplot(5,1,2);
%thresh_axes = gca;
%set(gca, 'Tag', 'threshaxes');
zscoreA = (g - mean(g))./std(g);  % in case Matlab zscore is not available
plot(abs(zscoreA));
ylabel('Std away from mean');
xlabel('Scan Number  -  horizontal axis for all plots');

thresh_x = 1:nscans;
thresh_y = z_thresh*ones(1,nscans);
line(thresh_x, thresh_y, 'Color', 'r');

%  Mark global intensity outlier images with vertical lines
axes_lim = get(gca, 'YLim');
axes_height = [axes_lim(1) axes_lim(2)];
for i = 1:length(glout_idx)
    line((glout_idx(i)*ones(1, 2)), axes_height, 'Color', 'r');
end

if realignfile == 1
	subplot(5,1,3);
    xa = [ 1:nscans];
	plot(xa,mv_data(:,1),'b-',xa,mv_data(:,2),'g-',xa,mv_data(:,3),'r-',...
       xa,mv_data(:,4),'r--',xa,mv_data(:,5),'b--',xa,mv_data(:,6),'g--');
    %plot(,'--');
	ylabel('ReAlignment');
	xlabel('Translation (mm) solid lines, Rotation (deg) dashed lines');
	legend('x mvmt', 'y mvmt', 'z mvmt','pitch','roll','yaw',0);
	h = gca;
	set(h,'Ygrid','on');
elseif realignfile == 0
    subplot(5,1,3);
	plot(mv0data(:,1:3));
	ylabel('Alignment (voxels)');
	xlabel('Scans. VERY APPROXIMATE EARLY-LOOK translation in voxels.');
	legend('x mvmt', 'y mvmt', 'z mvmt',0);
	h = gca;
	set(h,'Ygrid','on');
end 

subplot(5,1,4);   % Rapid movement plot
plot(delta);
ylabel('Motion (mm/TR)');
xlabel('Scan to scan movement (~mm). Rotation assumes 65 mm from origin');
y_lim = get(gca, 'YLim');
legend('Fast motion',0);
h = gca;
set(h,'Ygrid','on');

thresh_x = 1:nscans;
thresh_y = mv_thresh*ones(1,nscans);
line(thresh_x, thresh_y, 'Color', 'r');
   
% Mark all movement outliers with vertical lines
subplot(5,1,4)
axes_lim = get(gca, 'YLim');
axes_height = [axes_lim(1) axes_lim(2)];
for i = 1:length(mvout_idx)
    line((mvout_idx(i)*ones(1,2)), axes_height, 'Color', 'r');
end

%keyboard;
h_rangetext = uicontrol(gcf, 'Units', 'characters', 'Position', [10 10 18 2],...
        'String', 'StdDev of data is: ', 'Style', 'text', ...
        'HorizontalAlignment', 'left', ...
        'BackgroundColor', [0.8 0.8 0.8]);
h_rangenum = uicontrol(gcf, 'Units', 'characters', 'Position', [29 10 10 2], ...
        'String', num2str(gsigma), 'Style', 'text', ...
        'HorizontalAlignment', 'left',...
        'Tag', 'rangenum',...
        'BackgroundColor', [0.8 0.8 0.8]);
h_threshtext = uicontrol(gcf, 'Units', 'characters', 'Position', [25 8 16 2],...
        'String', 'Current threshold (std devs):', 'Style', 'text', ...
        'HorizontalAlignment', 'center', ...
        'BackgroundColor', [0.8 0.8 0.8]);
h_threshnum = uicontrol(gcf, 'Units', 'characters', 'Position', [44 8 10 2],...
        'String', num2str(z_thresh), 'Style', 'text', ...
        'HorizontalAlignment', 'left', ...
        'BackgroundColor', [0.8 0.8 0.8],...
        'Tag', 'threshnum');
h_threshmvtext = uicontrol(gcf, 'Units', 'characters', 'Position', [106 8 18 2],...
        'String', 'Motion threshold  (mm / TR):', 'Style', 'text', ...
        'HorizontalAlignment', 'center', ...
        'BackgroundColor', [0.8 0.8 0.8]);
h_threshnummv = uicontrol(gcf, 'Units', 'characters', 'Position', [126 8 10 2],...
        'String', num2str(mv_thresh), 'Style', 'text', ...
        'HorizontalAlignment', 'left', ...
        'BackgroundColor', [0.8 0.8 0.8],...
        'Tag', 'threshnummv');
h_threshtextpct = uicontrol(gcf, 'Units', 'characters', 'Position', [66 8 16 2],...
        'String', 'Current threshold (% of mean):', 'Style', 'text', ...
        'HorizontalAlignment', 'center', ...
        'BackgroundColor', [0.8 0.8 0.8]);
h_threshnumpct = uicontrol(gcf, 'Units', 'characters', 'Position', [86 8 10 2],...
        'String', num2str(z_thresh*pctmap), 'Style', 'text', ...
        'HorizontalAlignment', 'left', ...
        'BackgroundColor', [0.8 0.8 0.8],...
        'Tag', 'threshnumpct');
h_deweightlist = uicontrol(gcf, 'Units', 'characters', 'Position', [150 6 1 1 ],...
        'String', int2str(outdw_idx), 'Style', 'text', ...
        'HorizontalAlignment', 'left', ...
        'BackgroundColor', [0.8 0.8 0.8],...
        'Tag', 'deweightlist');
h_clipmvmtlist = uicontrol(gcf, 'Units', 'characters', 'Position', [152 6 1 1 ],...
        'String', int2str(clipout_idx), 'Style', 'text', ...
        'HorizontalAlignment', 'left', ...
        'BackgroundColor', [0.8 0.8 0.8],...
        'Tag', 'clipmvmtlist');
h_indextext = uicontrol(gcf, 'Units', 'characters', 'Position', [10 3 15 2],...
        'String', 'Outlier indices: ', 'Style', 'text', ...
        'HorizontalAlignment', 'left', ...
        'BackgroundColor', [0.8 0.8 0.8], ...
        'ForegroundColor', 'r');
h_indexedit = uicontrol(gcf, 'Units', 'characters', 'Position', [25 3.25 40 2],...
        'String', int2str(out_idx), 'Style', 'edit', ...
        'HorizontalAlignment', 'left', ...
        'Callback', 'art_outlieredit',...
        'BackgroundColor', [0.8 0.8 0.8],...
        'Tag', 'indexedit');
h_indexinst = uicontrol(gcf, 'Units', 'characters', 'Position', [66 3 40 2],...
        'String', '[Hit return to update after editing]', 'Style', 'text',...
        'HorizontalAlignment', 'left', ...
        'BackgroundColor', [0.8 0.8 0.8]);
h_deweighttext = uicontrol(gcf, 'Units', 'characters', 'Position', [115 1 21 2],...
        'String', 'Click to add margins for deweighting', 'Style', 'text', ...
        'HorizontalAlignment', 'center', ...
        'BackgroundColor', [0.8 0.8 0.8]);
h_clipmvmttext = uicontrol(gcf, 'Units', 'characters', 'Position', [94 1 21 2],...
        'String', 'Mark > 3mm movments for repair', 'Style', 'text', ...
        'HorizontalAlignment', 'center', ...
        'BackgroundColor', [0.8 0.8 0.8]);
if realignfile == 1
   h_repairtext = uicontrol(gcf, 'Units', 'characters', 'Position', [137 1 17 2],...
        'String', 'Writes repaired volumes', 'Style', 'text', ...
        'HorizontalAlignment', 'center', ...
        'BackgroundColor', [0.8 0.8 0.8]);
else  % realignfile == 0
   h_repairtext = uicontrol(gcf, 'Units', 'characters', 'Position', [137 1 17 2],...
        'String', 'WARNING!! DATA NOT REALIGNED', 'Style', 'text', ...
        'HorizontalAlignment', 'center', ...
        'BackgroundColor', [0.8 0.8 0.8]);
end
h_addmargin = uicontrol(gcf, 'Units', 'characters', 'Position', [120 3.25 10 2],...
        'String', 'Margin', 'Style', 'pushbutton', ...
        'Tooltipstring', 'Adds margins to deweight in estimation', ...
        'BackgroundColor', [ 0.7 0.9 0.7], 'ForegroundColor','k',...
        'Callback', 'art_addmargin');
h_clipmvmt = uicontrol(gcf, 'Units', 'characters', 'Position', [100 3.25 10 2],...
        'String', 'Clip', 'Style', 'pushbutton', ...
        'Tooltipstring', 'Marks displacements > 3 mm for repair', ...
        'BackgroundColor', [ 0.7 0.7 0.8 ], 'ForegroundColor','k',...
        'Callback', 'art_clipmvmt');
h_repair = uicontrol(gcf, 'Units', 'characters', 'Position', [140 3.25 10 2],...
        'String', 'REPAIR', 'Style', 'pushbutton', ...
        'Tooltipstring', 'Writes repaired images', ...
        'BackgroundColor', [ 0.9 0.7 0.7], 'ForegroundColor','r',...
        'Callback', 'art_repairvol');
h_up = uicontrol(gcf, 'Units', 'characters', 'Position', [10 8 10 2],...
        'String', 'Up', 'Style', 'pushbutton', ...
        'TooltipString', 'Raise threshold for outliers', ...
        'Callback', 'art_threshup');
h_down = uicontrol(gcf, 'Units', 'characters', 'Position', [10 6 10 2],...
        'String', 'Down', 'Style', 'pushbutton', ...
        'TooltipString', 'Lower threshold for outliers', ...
        'Callback', 'art_threshdown');

%guidata(gcf, g);
guidata(gcf,[g delta mv_data]);
setappdata(h_repair,'data',P);
setappdata(h_addmargin,'data2',repair1_flag);
setappdata(h_clipmvmt,'data3',repair1_flag);
%setappdata(h_repair,'data3',GoRepair);

% For GUI or RepairAlone, add deweighting margins to top plot
% Don't apply margins when motion adjustment was used.(v2.2)
if GoRepair == 0 | GoRepair == 1
    art_clipmvmt;
    art_addmargin;
end

if GoRepair == 1 | GoRepair == 2
    ImageFile = P(1,:);
    [imgpath, imagname] = fileparts(ImageFile); 
    [subjectpath, imagfold ] = fileparts(imgpath);
    [uppath2, up1name ] = fileparts(subjectpath);
    [uppath3, up2name ] = fileparts(uppath2);
    [uppath4, up3name ] = fileparts(uppath3);
    titlename = [ up3name, ' / ', up2name,' / ', up1name,' / ', imagfold ];
    gcf; subplot(5,1,1);
    title(titlename);
    %figname = ['artglobal', subname, '.jpg'];
    figname = ['artglobal_', up1name,'_', imagfold, '.jpg']; % added underscores to figure name (MemoLab)
    temp = pwd;
    cd(subjectpath);
    try
        saveas(gcf,figname);
    catch
        disp('Warning from art_global: Could not print figure');
    end
    cd(temp);
    % Return to directory in use before writing the jpg.
    art_repairvol(P);
end

if GoRepair == 4 
    
    % added same figure save even when no repairs are done (MemoLab)
    ImageFile = P(1,:);
    [imgpath, imagname] = fileparts(ImageFile); 
    [subjectpath, imagfold ] = fileparts(imgpath);
    [uppath2, up1name ] = fileparts(subjectpath);
    [uppath3, up2name ] = fileparts(uppath2);
    [uppath4, up3name ] = fileparts(uppath3);
    titlename = [ up3name, ' / ', up2name,' / ', up1name,' / ', imagfold ];
    gcf; subplot(5,1,1);
    title(titlename);
    %figname = ['artglobal', subname, '.jpg'];
    figname = ['artglobal_', up1name,'_', imagfold, '.jpg'];
    temp = pwd;
    cd(subjectpath);
    try
        saveas(gcf,figname);
    catch
        disp('Warning from art_global: Could not print figure');
    end
    cd(temp);
    
    %art_clipmvmt;
    outall_idx = unique([out_idx]);
    tempd = pwd;
    cd(fileparts(P(1,:)));
    save art_suspects.txt outall_idx -ascii
    %art_addmargin;
    cd(tempd)
end

% Define output variables (MemoLab)
delta_run    = delta;                % frame displacement at each TimePoint
art_suspects = unique(out_idx);      % list of all bad timepoints
FD_run       = mean(delta);          % mean frame displacement







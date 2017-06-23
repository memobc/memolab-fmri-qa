function [b, snrmean, snrmean_mid, spatialsnr_mean] = compute_snr(b, figFormat, pdfReport)
% COMPUTE_SNR
%
% This function generates SNR maps and generates and saves plots of
% intensity over time and mean SNR by run extended from scripts from here:
% 
% http://dbic.dartmouth.edu/wiki/index.php/Noise_Detection
%
% Author: Maureen Ritchey, Shao Fang Wang, Liang-Tien (Frank) Hsieh,
% July 2014
% Intended for use with memolab_batch_qa

%% CHECK FOR EXISTING FILES - MIGHT WANT TO FORCE IT TO RUN, remove option to turn off (buggy on next steps)
% Check existance of temporal files
fname1 = fullfile(b.dataDir, b.QAdir, [b.runs{1} '_SNR.nii']);
if exist(fname1, 'file');% check only for first run
    if b.auto_accept
        response_snr = 'y';
    else
        response_snr = input('SNR nii files exist. Do you want to overwrite the existing files? y/n \n','s');
    end
    if strcmp(response_snr,'y') == 1;
        disp('Continuing running compute_snr')
    else
        disp('Skipping compute_snr...')
        snrmean         = nan(length(b.runs));
        snrmean_mid     = nan(length(b.runs));
        spatialsnr_mean = nan(length(b.runs));
        return
    end
end

%% INITIALIZE VARIABLES ETC
snrmean = []; snrmean_mid = [];
close all
f  = figure('Visible','off');
st = figure('Visible','off');
g1 = figure('Visible','off');
g2 = figure('Visible','off');
g3 = figure('Visible','off');
s  = figure('Visible','off');
dv = figure('Visible','off');
gs = figure('Visible','off');
sp = figure('Visible','off');

%% LOOP OVER RUNS
% Start loop through runs
    
%Load brain mask
M    = spm_vol(b.mask);
mask = spm_read_vols(M);

%Initialize global signal suspects
GS_suspects = {};

%Loop over runs
for j = 1:length(b.runs)
    fprintf('\n%s (%02d / %02d runs)\n', b.runs{j}, j, length(b.runs))
    
    %% DATA HANDLING
    
    %Load preprocessed BOLDS
    disp('Preparing images')
    P     = b.rundir(j).rfiles;
    files = spm_vol(P);
    
    %Initialize variables
    outputname  = b.runs{j};% current run name
    data        = spm_read_vols(files(1));
    avg         = zeros(size(data));% create an empty variable sized by the image dimensions
    sd_tmp      = zeros(size(data));% again, creating empty variable sized by the image dimensions
    ntimepoints = size(files, 1);% figure out how many files you're dealing w
    slicematrix = [];
    
    %Store data and calculate metrics
    disp('Loading images')
    for i = 1:ntimepoints % loop across .nii files
        % concatenate all images into a 4 D matrix. The 4th dimension is time (n).
        data(:,:,:,i) = spm_read_vols(files(i));
    end
    avg = sum(data, 4)/ntimepoints; % average across time points
    
    %% TEMPORAL SNR
    disp('Computing temporal SNR')
    for i = 1:ntimepoints
        sd_tmp(:,:,:,i) = (avg-double(data(:,:,:,i))).^2;
        % average slices over within-brain voxels
        slicematrix = [ slicematrix squeeze(mean(mean(mask.*data(:,:,:,i))))];
    end
    sd_tmp = sum(sd_tmp,4);
    
    sd             = sqrt(sd_tmp/(ntimepoints-1));% standard SD calculation
    snr            = avg./sd;% standard SNR calculation
    meantimecourse = mean(slicematrix);
    
    %Clean up SNR variable
    snr             = mask.*snr; % mask out non-brain data
    snr(snr>2000)   = NaN; % mask out abnormal values
    snr(snr==0)     = NaN;
    snrmean         = [snrmean nanmean(nanmean(nanmean(snr(:,:,:))))];
    snrmean_mid     = [snrmean_mid nanmean(nanmean(snr(:,:,(ceil(size(data,3)/2)))))];
    snr(isnan(snr)) = 0;
    
    %output files to the QA dir
    outputdir  = fullfile(b.dataDir,b.QAdir);
    avg_output = files(1);
    sd_output  = files(1);
    snr_output = files(1);
    avg_output.fname    = fullfile(outputdir, [outputname '_average.nii']);
    sd_output.fname     = fullfile(outputdir, [outputname '_SD.nii']);
    snr_output.fname    = fullfile(outputdir, [outputname '_SNR.nii']);
    %     spm_write_vol(avg_output, avg); % typically don't need these but can uncomment if you want them
    %     spm_write_vol(sd_output, sd);
    spm_write_vol(snr_output, snr);
    
    %% SPATIAL SNR
    % Compute spatial SNR
    fprintf('Computing spatial SNR\n')
    clear signal_mean noise_SD;
    tmpdata = reshape(data,numel(data(:,:,:,1)),size(data,4))';
    tmpmask = reshape(mask,numel(data(:,:,:,1)),1)';
    signaldata = tmpdata; noisedata = tmpdata;
    signaldata(:,tmpmask==0 | sum(tmpdata)==0) = []; %remove voxels outside the brain
    noisedata(:,tmpmask==1 | sum(tmpdata)==0)  = []; %keep voxels outside the brain
    signal_mean = mean(signaldata'); %will return a row vector w/ the number of elements == number of images
    noise_SD    = std(double(noisedata')); %change to a double variable type (for increased precision) and take SD
    spatial_SNR = 0.655*(signal_mean./noise_SD);%constant (from the web) for fmri. then divide every element of the vector by SD_noise
    spatial_SNR_runs{j} = spatial_SNR;
    clear tmpdata tmpmask signaldata noisedata
    
    %% MAKE PLOTS
    %Save jpgs of slices
    avg_slices = save_slices(avg, fullfile(outputdir, [outputname '_average.jpg']), 0);
    sd_slices  = save_slices(sd, fullfile(outputdir, [outputname '_SD.jpg']), 0);
    snr_slices = save_slices(snr, fullfile(outputdir, [outputname '_tSNR.jpg']), 0);
    
    %plot slice x timepoint matrix
    figure(st)
    subplot(length(b.runs), 1, j);
    imagesc(slicematrix)
    xlabel('timepoints'); ylabel('slices');
    title(['Slice x Timepoint: ', b.curSubj, ' ', b.runs{j}], 'Interpreter', 'none');
    colorbar;
    
    %plot average slices
    figure(g1)
    subplot(length(b.runs), 1, j);
    image(avg_slices)
    colormap(gray)
    axis off
    title(['Mean: ' b.curSubj ' ' b.runs{j}], 'Interpreter', 'none');
    
    %plot sd slices
    figure(g2)
    subplot(length(b.runs), 1, j);
    image(sd_slices)
    colormap(jet)
    axis off
    title(['SD: ' b.curSubj ' ' b.runs{j}], 'Interpreter', 'none');
    
    %plot snr slices
    figure(g3)
    subplot(length(b.runs), 1, j);
    image(snr_slices)
    colormap(gray)
    axis off
    title(['tSNR: ', b.curSubj, ' ', b.runs{j}], 'Interpreter', 'none');

    clear data avg snr sd files
end

%% SAVE PLOTS

%finish timecourse plots
disp('Saving plots')
save_figs(st, 'timecourse_slice_plot', b, figFormat, pdfReport);
save_figs(g1, 'slices_avg_plot',       b, figFormat, pdfReport);
save_figs(g2, 'slices_sd_plot',        b, figFormat, pdfReport);
save_figs(g3, 'slices_tsnr_plot',      b, figFormat, pdfReport);

%plot mean SNR by run
figure(s)
snrmean = [snrmean' snrmean_mid'];
bar(snrmean);
title([b.curSubj '-Temporal SNR'], 'FontSize', 14);
xlabel('run'); ylabel('mean SNR');
legend('whole-brain', 'midslice', 'Location', 'EastOutside');
save_figs(s, 'tsnr_plot', b, figFormat, pdfReport, 6, 6); % square plot

% % Plot spatial SNR
figure(sp)
maxrunlen = size(spatial_SNR, 2);
%colors = {'r','g','b','k','y','m','c'};
colors = num2cell(hsv(size(spatial_SNR_runs,2)), 2);
for j = 1:size(spatial_SNR_runs,2) % coded this way in case runs have different lengths
    if size(spatial_SNR_runs{j}, 2) > maxrunlen
        maxrunlen = size(spatial_SNR_runs{j}, 2);
    end
    plot(spatial_SNR_runs{j}, 'Color', colors{(mod(j-1,length(colors))+1)}); % cycle through colors
    hold on
end
hold off
legend(b.runs, 'Location', 'EastOutside')
title([b.curSubj '-Spatial SNR'], 'FontSize', 14);
ylabel('Spatial SNR', 'FontSize', 14)
xlabel('Timepoint', 'FontSize', 14)
axis([0 maxrunlen -2 2]); axis 'auto y';
save_figs(sp, 'spatialSNR', b, figFormat, pdfReport, 6, 6); % square plot

%Compute mean spatial_SNR for summary info
for j = 1:size(spatial_SNR_runs, 2)
    spatialsnr_mean(j) = mean(spatial_SNR_runs{j}); % mean(spatial_SNR_runs,2);
end


end %snr_temporal

%% SUBFUNCTIONS
function [allslices] = save_slices(data, outputFile, writeFlag)
%this function takes a 3-D image and horizontally concatenates the select
%slices so they can be viewed in a row; returns the data for use in a plot,
%and optionally saves the multi-slice image
allslices = [];
for j = 1:3:size(data,3) % shows every third slice
    slice     = data(:,:,j)';
    allmin    = nanmin(nanmin(nanmin(data(:,:,:))));
    allmax    = nanmax(nanmax(nanmax(data(:,:,:))));
    slice     = round(63*(slice-allmin)/allmax)+1;
    allslices = [allslices slice];
end
if writeFlag,imwrite(allslices, gray, outputFile);end

end %save_slices

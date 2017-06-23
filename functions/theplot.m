function theplot(varargin)
% THEPLOT Create Power (2016)'s "The Plot" aka voxplot
%
%   Designed to work with Memolab QA Routine.
%   Adapted from Power (2016)'s demo.
%   Written by Kyle Kurkela, kyleakurkela@gmail.com August, 2016
%
%   Useage:
%
%   theplot(seg_image, bold_image, outdir, fd);
%
%   Plots the framewise displacement vector 'fd' above all of the voxels
%   contained in the timeseries 'bold_image', arranged by timepoint and
%   the first three compartments of spm12's segementation routine gray
%   matter, white matter, cerebrospinal fluid (i.e., c1*, c2*, and c3*
%   images). Where:
%
%   seg_image  = full path to an already segmented anatomical or mean
%                functional image. Function assumes the segmented images
%                are contained within the same directory as this image.
%
%   bold_image = full path to a 4-D timeseries .nii file that will be
%                displayed in the plot.
%
%   outdir     = full path to the directory where a jpg of theplot will be
%                saved.
%
%   fd         = a vector of framewise displacement values that correspond
%                to each volume of the bold_image nii timeseries.
%
%   For more information, see:
%
%   Power (2016) doi: dx.doi.org/10.1016/jneuroimage.2016.08.009
%   <a href="matlab:
%   web('http://www.jonathanpower.net/2016-ni-the-plot.html')">Jonathan Power's The Plot Demo</a>
%
%   See also memolab_batch_qa, run_theplot, spm_preproc

%====================================================================================
%			Step 0: Expand input arguments
%====================================================================================

% Expand Input Arguments
seg_image  = varargin{1}; % fullpath to this subject's segmented image
bold_image = varargin{2}; % fullpath to the this run's 4-D bold timeseries
outdir     = varargin{3}; % fullpath to the desired output directory
fd         = varargin{4}; % a vector of the framewise displacement measurments for each image of this timeseries
run        = varargin{5}; % a string which represents the name of the current functional run

[segDir, ~, ~]  = fileparts(seg_image);

%====================================================================================
%			Step 1: Can this script be run?
%====================================================================================
% We need to determine if the script can run. This script needs:
% The spm_segment

% Does segDir exist?
if ~exist(segDir, 'dir')
    error('The submitted segmented image directory does not seem to exist')
end

% Do the segmented c1*, c2*, and c3* images exist?
if isempty(spm_select('FPList', segDir, '^c1')) || isempty(spm_select('FPList', segDir, '^c2')) || isempty(spm_select('FPList', segDir, '^c3'))
    error('One of the segmented images c1*, c2*, c3* does not exist in this data directory. Segmentation may need to be run')
end

% Prepare the output directory
if ~exist(outdir, 'dir')
    mkdir(outdir)
end

% Does the fd vector match the length of the bold time series?
if size(bold_image,1) ~= length(fd)
    fprintf('\nImage Frames:\n\n')
    disp(bold_image)
    fprintf('\nLength of the FD vector: %0.0f\n\n', length(fd))
    assert(size(bold_image,1) == length(fd), 'The number of bold image frames does NOT match the number of framewise displacement values. See above.')
end

%====================================================================================
%			Step 2: Load, mask, and manipluate various images
%====================================================================================
% Loading masks, applying masks, and other image manipulation

%-- load the BOLD timeseries in RAS orientation

fname.bold = bold_image;
V          = spm_vol(fname.bold);
img.bold   = spm_read_vols(V);

%-- load the compartment masks in RAS orientation

% Grey Matter
fname.gm     = spm_select('FPList', segDir, '^c1');
V            = spm_vol(fname.gm);
mask.gm      = spm_read_vols(V);
mask.gm      = ~~mask.gm;

% White Matter
fname.wm     = spm_select('FPList', segDir, '^c2');
V            = spm_vol(fname.wm);
mask.wm      = spm_read_vols(V);
mask.wm      = ~~mask.wm;

% CSF
fname.csf   = spm_select('FPList', segDir, '^c3');
V           = spm_vol(fname.csf);
mask.csf    = spm_read_vols(V);
mask.csf    = ~~mask.csf;

%-- Misc Image Manipulations

% flatten and detrend the bold image
dd = size(img.bold);
im = reshape(img.bold, [dd(1)*dd(2)*dd(3) dd(4)]);

% mean and trend regressors
r0 = ones(1,dd(4));
r1 = linspace(0,1,dd(4));

% remove those terms
im = myregress([r0;r1], im);

%--Create the basic plot

img.masterflat = im; % the entire BOLD time series, entire brain flattened

gm      = img.masterflat(mask.gm(:),:);     % first grey matter voxels flattened
wm      = img.masterflat(mask.wm(:),:);     % first white matter voxels flattened
csf     = img.masterflat(mask.csf(:),:);    % first CSF voxels flattened

sz.gm   = size(gm,1);        % number of voxels in first grey matter
sz.wm   = size(wm,1);        % number of voxels in first white matter
sz.csf  = size(csf,1);       % number of voxels in first CSF

redlimz = sz.gm; % the total number of grey matter voxels

allmat  = [gm; wm; csf;]; % matrix of all voxels


%====================================================================================
%			Step 4: Plot the Plot
%====================================================================================
% Actually Ploting The Plot

%-- Set Up MATLAB figure

close all;
h = figure;
set(h, 'position', [100 100 1920 1080]);
set(h, 'visible', 'on');

%-- Figure settings (e.g., colors)

sett.gmclr  = [0 .5 1];
sett.wmclr  = [0 1 0];
sett.csfclr = [1 1 0];
gmlimz      = [-20 20];
xlimz       = [1 dd(4)];


%-- Head Motion (fd)

subplot(5,1,1);
plot(1:dd(4), fd, 'r');
xlim([1 dd(4)]);
if max(fd) > 1
    ylim([0 max(fd)]);
    set(gca, 'ytick', 0:1:max(fd));
else
    ylim([0 1]);
    set(gca, 'ytick', 0:.5:1);
end
title('Head motion');
ylabel('Framewise Displacement (mm)');


%--The Heatmap

subplot(5, 1, 2:5);

imagesc(allmat, gmlimz);
set(gca, 'ytick', []);
colormap(gray);
hh = hline(redlimz(1), 'g');
set(hh, 'linewidth' , 2);
xlim(xlimz);
title('Timeseries of in-brain voxels, arranged by SPM NewSegment generated Masks');
ylabel('Voxels');
xlabel('Volume #');

%--This side color bands

% First Grey Matter Compartment
stx = 0; endx = 2;
sty = 0; endy = sz.gm;
v = [stx sty; stx endy; endx endy; endx sty];
f = [1 2 3 4];
patch('Faces', f, 'Vertices', v, 'FaceColor', .4*sett.gmclr);

% First White Matter Compartment
sty = sz.gm; endy=sz.gm+sz.wm;
v = [stx sty; stx endy; endx endy; endx sty];
f = [1 2 3 4];
patch('Faces', f, 'Vertices', v, 'FaceColor', .8*sett.wmclr);

% Cerbro-Spinal Fluid Compartment
sty = sz.gm+sz.wm; endy=sz.gm+sz.wm+sz.csf;
v = [stx sty; stx endy; endx endy; endx sty];
f =[ 1 2 3 4];
patch('Faces', f, 'Vertices', v, 'FaceColor', sett.csfclr);


%--Far side color bands

% First Grey Matter Compartment
stx = xlimz(end)-1; endx=xlimz(end)+1;
sty = 0; endy=sz.gm;
v = [stx sty; stx endy; endx endy; endx sty];
f = [1 2 3 4];
patch('Faces', f, 'Vertices', v, 'FaceColor', .4*sett.gmclr);

% First White Matter Compartment
sty = sz.gm; endy=sz.gm+sz.wm;
v = [stx sty; stx endy; endx endy; endx sty];
f = [1 2 3 4];
patch('Faces', f, 'Vertices', v, 'FaceColor', .8*sett.wmclr);

% Cerbro-Spinal Fluid Compartment
sty = sz.gm+sz.wm; endy=sz.gm+sz.wm+sz.csf;
v = [stx sty; stx endy; endx endy; endx sty];
f = [1 2 3 4];
patch('Faces', f, 'Vertices', v, 'FaceColor', sett.csfclr);


%--Save plot

oname = fullfile(outdir, [run '_voxplot.png']);
set(h, 'paperpositionmode', 'auto');
print(gcf, '-dpng', oname);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [resid, pred, b]= myregress(r, tc, varargin)
        
        % presumes vox x time input variable structure
        
        if isempty(varargin)
            % use all timepoints
            b     = r'\tc';
            pred  = r'*b;
            resid = tc-pred';
            
        else
            % use only specified timepoints
            tmask = varargin{1,1};
            b     = r(:,tmask)'\tc(:,tmask)';
            pred  = r'*b;
            resid = tc-pred';
        end
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

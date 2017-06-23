function [] = memolab_batch_qa()
% The MemoLab's Batch QA script, which builds on prior DML QA scripts
%
% See README for full description of features and output. 
%
% Authors: Kyle Kurkela (BC), Maureen Ritchey (BC, Davis)
% Contributors: Shao-Fang Wang (Davis), Liang-Tien Hsieh (Davis),
% and Halle Zucker (Davis), August 2016

%====================================================================================
%			Specify Variables
%====================================================================================

%-- File Type
% Are you using DCM or NII files? If NII, they should be unprocessed
% .*\.nii files. Alternatively they can be a 4D .nii file (single file per
% run; can be .*\.nii or .*\.nii.gz, although gz files will be unzipped for
% SPM). 'DCM' or 'NII'. Note: this QA routine is NOT compatible with
% .img/.hdr. Please convert .img/.hdr to .nii prior to running routine.

fileType    = 'NII';

%-- Directory Information
% Paths to relevant directories.
% dataDir   = path to the directory that houses the MRI data
% scriptdir = path to directory housing this script (and auxiliary scripts)
% QAdir     = Name of output QA directory

dataDir     = '/fullpath/to/the/data';
scriptdir   = '/fullpath/to/this/script'; % fileparts(mfilename('fullpath'));
QAdir       = 'QA';

%-- Info for Subjects
% Subject-specific information.
% subjects  = cellstring containing the subject IDs
% runs      = cellstring containg the IDs for each BOLD time series
%
% Assumes that all files have unique filenames that can be identified with
% a combination of the cell strings above. For example, bold files NEED to
% look something like:
%   /dataDir/sub-001/func/sub-001_encoding_run-001_bold.nii
%   /dataDir/sub-001/func/sub-001_encoding_run-002_bold.nii
%   /dataDir/sub-001/func/sub-001_retrieval_run-001_bold.nii
%   /dataDir/sub-001/func/sub-001_retrieval_run-002_bold.nii
%
%  See BIDS format

subjects    = {'sub-001'};
runs        = {'encoding_run-01' 'encoding_run-02' ...
               'retrieval_run01' 'retrieval_run02'};

%-- Figure Format
% Format figures should be saved in. Options: '-pdf' or '-png'
% pdf requires ghostscript (http://pages.uoregon.edu/koch/; see export_fig documentation)
% but looks better.

figFormat   = '-png';

%-- PDF Report
% If pdf was selected above, you can choose to output only a single PDF report
% (1) rather than individual figures (0, default).

pdfReport   = 0;

%-- Frame Displacement Threshold for ArtRepair
% Power et al. 2012 uses .5 but 2014 paper recommends an even lower
% threshold if dataset permits it.

mv_thresh    = .2;

%-- Global Signal Intensity Threshold
% In terms of percent deviation from mean. ARTrepair default was 1.3.

Percent_thresh    = 1.3;

%-- Art Global Options
% HeadMaskType is a flag that tells the art global routine in Art Repair which brain mask to use.
%	1 = use SPM mask
% 	4 = use Automask - creates ArtifactMask.nii file
% RepairType is a flag that controls the functionality of the art global routine:
%	1 = for ArtifactRepair alone (0.5 movement and add margin)
%	2 = for Movement Adjusted images  (0.5 movement, no margin)
%	0 = No repairs are done, bad scans are found. Listed in art_suspects.txt for motion adjustment.

HeadMaskType = 4;
RepairType 	 = 0;

%-- Prefix for Spike Regressor File.
% Note: can be empty, i.e. ''

spikePrefix = '';

%-- Auto-accept
% Do you want to run all the way through without asking for user input?
% if 0: will prompt you to take action;
% if 1: skips realignment and ArtRepair if already run, overwrites output files

auto_accept = 0;

%-- Run The Plot Flag
% Do you want to create a voxplot for each run for each subject?
% if 0: will skip creating the plot
% if 1: will create the plot and segment the mean functional (if necessary)
% - this is meant to be a rough estimate for visualization purposes only

runtheplot = 1;

%====================================================================================
%			Routine (DO NOT EDIT BELOW THIS BOX!!)
%====================================================================================

%-- Clean up
close all
clc
fprintf('Initializing and checking paths.\n')

%-- Add paths
addpath(genpath(fullfile(scriptdir, 'functions')));
if exist(fullfile(scriptdir, 'vendor'),'dir')
    addpath(genpath(fullfile(scriptdir, 'vendor')));
end

%-- Check for required functions

% SPM
if exist('spm','file') == 0
    error('SPM must be on the path.')
end

% export_fig
if exist('export_fig','file') == 0
    error('export_fig must be on the path.')
end

% nan stats tools
if exist('nanmean','file') == 0
    error('nan tools must be on the path.')
end

% Art Repair
if exist('art_global','file') == 0
    error('ArtRepair must be on the path.')
end

if runtheplot % voxplot specific requirements
    % hline/vline
    if exist('vline', 'file') == 0
        error('hline and vline must be on the path.')
    end
end

fprintf('Running QA script. Results will be saved in QA directory: %s\n', QAdir)

%--Loop over subjects
for i = 1:length(subjects)
    
    % Define variables for individual subjects - General
    b.curSubj   = subjects{i};
    b.runs      = runs;
    b.dataDir   = fullfile(dataDir, b.curSubj);
    
    %%% Alternatively, if there is an initializeVars script set up, call that
    %%% see https://github.com/ritcheym/fmri_misc/tree/master/batch_system    
    %     b = initializeVars(subjects,i);
    
    % Define variables for individual subjects - QA General
    b.scriptdir   = scriptdir;
    b.QAdir       = QAdir;
    b.auto_accept = auto_accept;
    b.messages    = sprintf('Messages for subject %s:\n', subjects{i});
    
    % Define variables for individual subjects - Art Repair Specific
    b.mv          = mv_thresh;
    b.percent     = Percent_thresh;
    b.headMask    = HeadMaskType;
    b.repairType  = RepairType;
    
    fprintf('Subject %s \n', b.curSubj)
    
    % Check whether QA has already been run for a subject
    if ~exist(fullfile(b.dataDir, b.QAdir), 'dir')
        fprintf('Creating QA directory for subject %s\n', b.curSubj)
        mkdir(fullfile(b.dataDir, b.QAdir));
    else
        fprintf('QA directory exists for subject %s\n', b.curSubj)
        if b.auto_accept
            response_QA = 'y';
        else
            response_QA = input('Do you want to continue anyway? y/n \n', 's');
        end
        if strcmp(response_QA, 'y') == 1
            disp('Continuing running QA script')
            if pdfReport == 1
                % check if there's already a QA_report pdf file
                if exist(fullfile(b.dataDir, b.QAdir, ['QA_report_' b.curSubj '.pdf']),'file')
                    if b.auto_accept
                        response_QApdf = 'y';
                    else
                        response_QApdf = input('Do you want to discard the previous PDF report? y/n \n','s');
                    end
                    if strcmp(response_QApdf,'y') == 1
                        delete(fullfile(b.dataDir, b.QAdir, ['/QA_report_'  b.curSubj '.pdf']));
                    end
                end
            end
        else
            error('Skipping subject %s\n', b.curSubj)
        end
    end % if for checking QA directories
    
    % Initialize diary for saving output
    diaryname = fullfile(b.dataDir, b.QAdir, 'QA_diary_output.txt');
    diary(diaryname);
    
    % Convert dicom images or find nifti images
    if strcmp(fileType, 'DCM')
        fprintf('--Converting DCM files to NII format--\n')
        [b] = convert_dicom(b);
        fprintf('------------------------------------------------------------\n')
        fprintf('\n')
    elseif strcmp(fileType, 'NII')
        fprintf('--Finding NII files (no DCM conversion required)--\n')
        [b] = find_nii(b);
        fprintf('------------------------------------------------------------\n')
        fprintf('\n')
    else
        error('Specify a valid fileType (''DCM'' or ''NII'')');
    end
    
    % Run realignment to generate motion parameters
    fprintf('--Realigning images using spm_realign--\n')
    [b] = batch_spm_realign(b);
    fprintf('------------------------------------------------------------\n')
    fprintf('\n')
    
    % Run art_global
    fprintf('--Running Art Global--\n')
    [all_suspects, FD_mean_run, FD_run] = run_art_global(b);
    fprintf('------------------------------------------------------------\n')
    fprintf('\n')
   
    % Calculate temporal SNR & plot intensity changes
    fprintf('--SNR & Intensity Calculations--\n')
    art_automask(b.meanfunc, -1, 1);
    b.mask = spm_select('ExtFPListRec', b.dataDir, 'ArtifactMask.nii', 1);
    [b, snrmean, snrmean_mid, spatialsnr_mean] = compute_snr(b, figFormat, pdfReport);
    fprintf('------------------------------------------------------------\n')
    fprintf('\n')

    if runtheplot
        
        % Create Power (2016)'s "The Plot"
        fprintf('--Creating the voxplot--\n')
        run_theplot(b, FD_run)
        fprintf('------------------------------------------------------------\n')
        fprintf('\n')
    
    end
    
    % Generate spike regressors
    fprintf('--Generate spike regressors for subject %s --\n', b.curSubj)
    [b, numBad] = create_spike_regs(b, all_suspects, spikePrefix);
    fprintf('------------------------------------------------------------\n')
    
    % Save out summary information
    fprintf('--Save summary info for subject %s --\n', b.curSubj)
    save_sum_info(b, FD_mean_run, FD_run, snrmean, snrmean_mid, spatialsnr_mean, numBad);
    fprintf('=============================================================\n')
    
end % isub

fprintf('Done QA script\n')
diary off

end % main function
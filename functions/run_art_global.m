function [all_suspects, FD_mean_run, FD_run] = run_art_global(b)
% RUN_ART_GLOBAL wrapper script for running Art Repair's art_global routine
% over several functional runs
%
% Intended for use with Memo Lab QA Routine
% Written by Kyle Kurkela and Maureen Ritchey, August 2016
%
%   Usage:
%
%   run_art_global(b)
%
%   Runs the Memolab's tweaked version of Art Repair's art_global routine
%   over several functional runs. Where:
%
%   b = memolab batch structure containing the fields:
%
%       b.meanfunc = fullpath string to the mean functional images created during
%                    SPM's realignment procedure
%
%       b.dataDir  = fullpath string to the directory where the functional MRI data
%                   is being stored
%
%       b.QAdir    = string, name of the directory where QA output is stored
%
%       b.rundir   = a 1 x n structure array, where n is the number of
%                    runs, with the fields:
%
%           b.rundir(n).files = cellstring of full paths to each frame of
%                               this run's nii timeseries
%
%           b.rundir(n).rp    = full path string to this run's motion
%                               parameters as estimated by SPM's
%                               realignment routine
%
%       b.runs        = cellstring with the name of the directories containing
%                       each functional run.
%
%       b.headMask    = art_global's mask flag, see art_global,
%                       memolab_batch_qa
%
%       b.repairType  = art_global's repair type flag, see art_global,
%                       memolab_batch_qa
%
%       b.percent     = art_global's global signal intensity threshold, see
%                       art_global, memolab_batch_qa
%
%       b.mv          = art_global's movement threshold, see art_global,
%                       memolab_batch_qa
%
%   See also memolab_batch_qa, art_global, art_global_qa

%====================================================================================
%			% Step 1: Initalize Cell Arrays
%====================================================================================

% Initalize FD cell arrays
all_suspects = cell(1,length(b.rundir));
FD_mean_run  = cell(1,length(b.rundir));
FD_run       = cell(1,length(b.rundir));

%====================================================================================
%			% Step 2: Loop Art Global over runs
%====================================================================================

% Check if ArtRepair was already run and if so, whether it should be run
% again
runflag = 1;
if size(spm_select('FPlistRec', b.dataDir, '^art'), 1) > 0
    if b.auto_accept
        response = 'n';
    else
        response = input('ArtRepair appears to have been run already. Do you want to run it again? y/n \n','s');
    end
    if strcmp(response,'y')==1
        disp('Continuing running ArtRepair')
    else
        disp('Skipping ArtRepair')
        runflag = 0;
        response = input('Would you like to try and load the FD vectors from summary_info.mat? y/n \n','s');
        if strcmp(response, 'y')
            summary_info = [];
            summary_info_file = fullfile(b.dataDir,b.QAdir,[b.curSubj '_QA_summary_information.mat']);
            load(summary_info_file)
            FD_run = horzcat(summary_info{3,2:end});
        else
        end
    end
end

% Run Art Global for each run
if runflag
    for r = 1:length(b.rundir)
        fprintf('%s (%d / %d runs)\n', b.runs{r}, r, length(b.rundir))
        [art_suspects, FD_mean_run{r}, FD_run{r}] = art_global_qa(b.rundir(r).files, b.rundir(r).rp, b.headMask, b.repairType, b.percent, b.mv);
        fprintf('------------------------------------------------------------\n')
        fprintf('\n')
        all_suspects{r} = art_suspects;
    end
end

%====================================================================================
%			% Step 3: Move Art Global Output
%====================================================================================

% Move Art Global Figures to QA dir
if runflag
    fprintf('--Moving Art Global figures...--\n')
    movefile(fullfile(b.dataDir, 'artglobal*.jpg'), fullfile(b.dataDir, b.QAdir))
end

end
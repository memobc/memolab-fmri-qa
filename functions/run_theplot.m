function run_theplot(b, FD_run)
% RUN_THEPLOT Wrapper function for running theplot.m over multiple runs
%
%   Intended for use with Memolab QA Routine.
%   Written by Kyle Kurkela and Maureen Ritchey, August 2016.
%
%   Usage:
%
%   run_theplot(b, FD_run)
%
%   Creates Power 2016's "The Plot" for each inputted functional run using
%   theplot.m. Where:
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
%           b.rundir(n).rfiles = full path to run n's realigned 
%                                and resliced bold_nii time series
%
%       b.runs    = cellstring with the name of the directories containing
%                   each functional run.
%
%   FD_run = a cell array containing vectors of the framewise displacement
%            values for each timepoint in all functional runs, where 
%            FD_run{1} = framewise displacement vector for run 1, FD_run{2}
%            = framewise displacement vector for run 2, and so on.
%
%   See also memolab_batch_qa, theplot, newsegment, spm_jobman, spm_preproc
    
%====================================================================================
%			% Step 1: Segmentation Check
%====================================================================================

    % Does segmentation need to be run?
    [segDir, ~, ~] = fileparts(b.meanfunc);    
    if isempty(spm_select('FPList', segDir, '^c1')) || ...
       isempty(spm_select('FPList', segDir, '^c2')) || ...
       isempty(spm_select('FPList', segDir, '^c3'))
        
        fprintf('--Segmentation has not yet been run. Running segmentation...--\n')
        newsegment(b.meanfunc);
        fprintf('------------------------------------------------------------\n')
        fprintf('\n')
    end

%====================================================================================
%			% Step 2: Create the Plot for each Functional Run
%====================================================================================
    
    % QA dir full path
    full_path_to_QAdir = fullfile(b.dataDir, b.QAdir);
    
    % For each run...
    for r = 1:length(b.rundir)
        fprintf('%s (%d / %d runs)\n', b.runs{r}, r, length(b.rundir))
        theplot(b.meanfunc, b.rundir(r).rfiles, full_path_to_QAdir, FD_run{r}, b.runs{r})
        fprintf('------------------------------------------------------------\n')
        fprintf('\n')
    end
    
end
function [b] = find_nii(b)
% FIND_NII finds nii files in the run directory and stores them in the 'b' 
%          structure.
%
%   Usage:
%
%	b = find_nii(b)
%   
%   input arguments:
%
%	b = memolab qa batch structure containing the fields:
%
%       b.runs      = cellstring with IDs for each functional time series
%
%       b.dataDir   = fullpath string to the directory where the functional MRI data
%                     is being stored
%
%       b.rundir    = a 1 x n structure array, where n is the number of
%                     runs
%
%   output arguments:
%
%       b = memolab qa batch structure containg the new fields:
%
%       b.rundir(n).files = a character array containg the full paths to
%                           the recently converted .*\.nii files
%
%
%   Intended for use with Memolab QA Routine.
%   Written by Maureen Ritchey, circa 2014
%
%   See also: gunzip, spm_select

for irun = 1:length(b.runs)
    
    % Select this run's full nii timeseries using a regular expression
    regularExp           = [ '^' b.curSubj '.*' b.runs{irun} '.*bold\.nii$'];
    b.rundir(irun).files = spm_select('ExtFPListRec', b.dataDir, regularExp, Inf);
    
    % Check if files are found
    if size(b.rundir(irun).files, 1) > 0
        fprintf('%0.0f nii files found.\n', size(b.rundir(irun).files, 1));
    else
        % Check for a .gz file that matches the regular expression
        gzfiletestExp = [ '^' b.curSubj '.*' b.runs{irun} '.*bold\.nii\.gz$'];
        gzfile        = spm_select('FPListRec', b.dataDir, gzfiletestExp);
        
        if ~isempty(gzfile)
            fprintf('Did not find nii files. Unzipping gz files for SPM:\n');
            disp(gzfile)
            gunzip(gzfile)
            b.rundir(irun).files = spm_select('ExtFPListRec', b.dataDir, regularExp, Inf);
        else
            error('No nii or matching gz files found.\n');
        end
    end
    
end

end

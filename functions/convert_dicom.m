function [b] = convert_dicom(b)
% CONVERT_DICOM  converts dicom files and stores paths in the b structure
%
%   Intended for use with memolab_batch_qa.
%   Written by Maureen Ritchey (BC, Davis), circa 2014
%   
%   Usage:
%
%	b = convert_dicom(b)
%   
%   input arguments:
%
%	b = memolab qa batch structure containing the fields:
%
%       b.runs      = cellstring with the name of the directories containing
%                     each functional run
%
%       b.dataDir   = fullpath string to the directory where the functional MRI data
%                     is being stored
%
%       b.rundir    = a 1 x n structure array, where n is the number of
%                     runs
%
%       b.scriptdir = fullpath string to the directory where the memolab
%                     scripts
%
%   output arguments:
%
%   b = memolab qa batch structure containg the new fields:
%
%       b.rundir(n).files = a character array containg the full paths to
%                           the recently converted .*\.nii files
%
%
% See also: find_nii, spm_dicom_headers, spm_dicom_convert

for irun = 1:length(b.runs)
    
    rundir   = fullfile(b.dataDir, b.runs{irun});
    dcmfiles = spm_select('FPList', rundir, '.*dcm');
    fprintf('In directory %s:\n%0.0f dcm files found...\n', rundir, size(dcmfiles,1));
    
    % Check whether there are already nifti files
    niifiles = spm_select('FPList', rundir, '.*\.nii');
    if size(niifiles, 1) > 0
        fprintf('There are already %0.0f nii files in this folder.\n', size(niifiles,1));
        response = input('Are you sure you want to use dicom files? y/n \n', 's');
        if strcmp(response,'y') == 1
            disp('Continuing running dicom conversion')
        else
            error('Change fileType to ''NII''')
        end
    end
    
    % Convert dicom images
    dcmhdr    = spm_dicom_headers(dcmfiles);
    cd(rundir);
    dcmoutput = spm_dicom_convert(dcmhdr, 'all', 'flat', 'nii');
    
    b.rundir(irun).files = cell2mat(dcmoutput.files);
    fprintf('%0.0f files converted to nii.\n', size(b.rundir(irun).files, 1));
    cd(b.scriptdir);
    
end

end
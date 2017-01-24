function [] = newsegment(img)
% NEWSEGMENT run SPM's new segementation routine
%
% Function designed to run a standard SPM segementation routine using the 
% new segementation approach.
%
%   Usage:
%
%   newsegment(img)
%
%   Runs SPM's new segmentation routine using SPM's job manager. Where:
%
%       img = full path string to image that is to be segmented
%
%   See also spm_jobman, spm_preproc


    matlabbatch{1}.spm.tools.preproc8.channel.vols     = {img};
    matlabbatch{1}.spm.tools.preproc8.channel.biasreg  = 0.0001;
    matlabbatch{1}.spm.tools.preproc8.channel.biasfwhm = 60;
    matlabbatch{1}.spm.tools.preproc8.channel.write    = [0 0];
    matlabbatch{1}.spm.tools.preproc8.tissue(1).tpm    = {fullfile(fileparts(which('spm')), 'tpm', 'TPM.nii,1')};
    matlabbatch{1}.spm.tools.preproc8.tissue(1).ngaus  = 2;
    matlabbatch{1}.spm.tools.preproc8.tissue(1).native = [1 0];
    matlabbatch{1}.spm.tools.preproc8.tissue(1).warped = [0 0];
    matlabbatch{1}.spm.tools.preproc8.tissue(2).tpm    = {fullfile(fileparts(which('spm')), 'tpm', 'TPM.nii,2')};
    matlabbatch{1}.spm.tools.preproc8.tissue(2).ngaus  = 2;
    matlabbatch{1}.spm.tools.preproc8.tissue(2).native = [1 0];
    matlabbatch{1}.spm.tools.preproc8.tissue(2).warped = [0 0];
    matlabbatch{1}.spm.tools.preproc8.tissue(3).tpm    = {fullfile(fileparts(which('spm')), 'tpm', 'TPM.nii,3')};
    matlabbatch{1}.spm.tools.preproc8.tissue(3).ngaus  = 2;
    matlabbatch{1}.spm.tools.preproc8.tissue(3).native = [1 0];
    matlabbatch{1}.spm.tools.preproc8.tissue(3).warped = [0 0];
    matlabbatch{1}.spm.tools.preproc8.tissue(4).tpm    = {fullfile(fileparts(which('spm')), 'tpm', 'TPM.nii,4')};
    matlabbatch{1}.spm.tools.preproc8.tissue(4).ngaus  = 3;
    matlabbatch{1}.spm.tools.preproc8.tissue(4).native = [0 0];
    matlabbatch{1}.spm.tools.preproc8.tissue(4).warped = [0 0];
    matlabbatch{1}.spm.tools.preproc8.tissue(5).tpm    = {fullfile(fileparts(which('spm')), 'tpm', 'TPM.nii,5')};
    matlabbatch{1}.spm.tools.preproc8.tissue(5).ngaus  = 4;
    matlabbatch{1}.spm.tools.preproc8.tissue(5).native = [0 0];
    matlabbatch{1}.spm.tools.preproc8.tissue(5).warped = [0 0];
    matlabbatch{1}.spm.tools.preproc8.tissue(6).tpm    = {fullfile(fileparts(which('spm')), 'tpm', 'TPM.nii,6')};
    matlabbatch{1}.spm.tools.preproc8.tissue(6).ngaus  = 2;
    matlabbatch{1}.spm.tools.preproc8.tissue(6).native = [0 0];
    matlabbatch{1}.spm.tools.preproc8.tissue(6).warped = [0 0];
    matlabbatch{1}.spm.tools.preproc8.warp.mrf    = 0;
    matlabbatch{1}.spm.tools.preproc8.warp.reg    = 4;
    matlabbatch{1}.spm.tools.preproc8.warp.affreg = 'mni';
    matlabbatch{1}.spm.tools.preproc8.warp.samp   = 3;
    matlabbatch{1}.spm.tools.preproc8.warp.write  = [0 0];
    spm_jobman('initcfg');
    spm_jobman('run', matlabbatch);
    
end
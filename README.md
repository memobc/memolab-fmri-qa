# MemoLab fMRI QA

SPM-based scripts for doing general QA on fMRI data, including identification of suspect timepoints with ArtRepair and generation of several diagnostic plots. QA can be run on several subjects at once.

# Requirements:
- [Matlab](http://www.mathworks.com/products/matlab/) - tested with R2015a
- [SPM](http://www.fil.ion.ucl.ac.uk/spm/software/spm12/) - tested with SPM12
- [ArtRepair toolbox](http://cibsr.stanford.edu/tools/human-brain-project/artrepair-software.html) - tested with v5b3
- MATLAB's [export_fig toolbox](http://www.mathworks.com/matlabcentral/fileexchange/23629-export-fig)
- [hline and vline](http://www.mathworks.com/matlabcentral/fileexchange/1039-hline-and-vline)
- MATLAB's [Statistics and Machine Learning Toolbox](http://www.mathworks.com/products/statistics/) (or just nanmean, nanmin, nanmax, nanstd)
  - Note: MATLAB's Statistics and Machine Learning Toolbox should be provided in the base version of MATLAB in ~R2006+

Note that any toolboxes placed in a folder called `vendor` within the repository will be automatically added to the Matlab path (and ignored by GitHub).

# Features:
- Automatic detection of necessary requirements
- Automatic detection of whether QA has already been run for a given subject
- Converts DICOM images or finds already-converted Nifti images
- Runs realignment to generate motion parameters
- Runs ART repair to identify bad timepoints
 - Calculate and plot global average signal over time
 - Calculate and plot standard deviation away from grand mean over time
 - Plot motion parameters
 - Calculate and plot framewise displacement (fast motion)
 - Flag suspect volumes based on framewise displacement and/or global signal
- Calculates temporal SNR & plots intensity changes
- Generates spike regressors based on ART repair suspects
- Generates voxplots (Power 2016 "the plot") for each run
- Saves summary information for each subject

# Usage:
1. Edit the top section of `memolab_batch_qa.m` to specify your input paths, variables, etc. All inputs are described in detail within the script.
2. Run the `memolab_batch_qa` script. It will run all steps detailed below.
3. If you need to run the script again, you can do so and it will skip over steps that have already been run (if you want it to).


# Processing steps:
1. Add SPM path
2. Begin loop over subjects
3. Check whether QA has already been run for a subject by checking the existence of the resulting files – *gives user the option to stop if it’s been run*
4. Convert DICOM images or find Nifti images by SPM functions: `convert_dicom.m` or `find_nii.m`
5. Run realignment to generate motion parameters: `batch_spm_realign.m`
  - Check whether motion parameters file exists for only the first run – *gives user the option to skip this step if it’s already been run*
  - Save motion parameters file for each run: rp*.txt
6. Run ArtRepair to identify bad timepoints: `run_art_global.m`
7. Calculate temporal SNR and spatial SNR: `compute_snr.m`
  - Check whether files exist: (6 in total) – *gives user the option to overwrite*
  - Apply the ArtifactMask.nii generated from ArtRepair to mask out non-brain data
  - Temporal SNR
  - Spatial SNR
8. Create Power et al. 2016's voxplot: `run_theplot.m`
9. Create spike regressors for late addition to a GLM: `create_spike_regs.m`
10. Save out summary information for each subject (cell array format): `save_sum_info.m`
  - Frame displacement
  - Temporal SNR
  - Temporal SNR midslice
  - Spatial SNR

# Output:
## Plots
-	`slices_avg_plot_[subject].png` – plot of average signal in image slices; should look like a brain!
-	`slices_sd_plot_[subject].png` – plot of signal SD in image slices; should be mostly blue
-	`slices_tsnr_plot_[subject].png` – plot of temporal SNR in image slices
-	`tsnr_plot_[subject].png` – plot of mean temporal SNR values for each run, both for the whole brain and the middle slice
-	`spatialSNR_[subject].png` – plot of spatial SNR over time
-	`timecourse_slice_plot_[subject].png` – plot of signal intensity over time by each slice; should look like a rainbow without any obviously deviant individual slices or timepoints
- `artglobal_[subject]_[run].png` - plot of the global mean signal, standard devation away from the mean signal, the 6 movement parameters and the framewise displacement for each timepoint. Bad timepoints are flagged with vertical lines
- `[run]_voxplot.png` - Power et al. 2016's "The Plot" for this run. Framewise displacement is plotted at the top; every voxel in a timepoints is represented by a single column in the body of the plot; the top 1/3 of the body if the grey matter (blue), then white matter (green), then CSF (yellow); there is a solid green line denoting the change from Grey Matter to White Matter

## Other files
-	`[subject]_QA_summary_information.mat` – contains subject summary info
-	`[run]_SNR.nii` – NII of temporal SNR
- `QA_diary_output.txt` - record of the MATLAB command window while running the memolab_MRI_qa pipeline
- `QA_report_[subject].pdf` - optional report that concatenates most of the plots, putting them into a single PDF (in batch script, set `pdfReport=1`)

# More information:
## References
- ArtRepair: https://www.nitrc.org/projects/art_repair/
- ThePlot aka voxplot: Power. NeuroImage (2016) doi: 10.1016/j.neuroimage.2016.08.009.

## Authors & Contributors
- DML QA scripts: Maureen Ritchey, Shao-Fang Wang, Liang-Tien Hsieh, and Halle R. Zucker, July 2014  
- MemoLab QA toolbox: Maureen Ritchey and Kyle Kurkela, August 2016

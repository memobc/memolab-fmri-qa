function [] = save_sum_info(b, FD_mean_run, FD_run, snrmean, snrmean_mid, spatialsnr_mean, numBad)
% SAVE_SUM_INFO  saves out summary information for memolab_batch_qa
%
%   Usage:
%
%	save_sum_info(b, FD_mean_run, FD_run, snrmean, snrmean_mid, spatialsnr_mean, numBad)
%   
%   input arguments:
%
%       b = memolab qa batch structure containing the fields:
%
%           b.runs      = cellstring with IDs for each functional time
%                         series
%
%           b.curSubj   = string, the ID for the current subject
%
%           b.dataDir   = fullpath string to the directory where the functional MRI data
%                         is being stored
%
%           b.QAdir     = string, cotaining the name of the output QA directory
%
%           b.messages  = various messages collected throughout the QA script
%
%
%       FD_mean_run     = a cell array containg the mean framewise displacement for
%                         each run
%
%       FD_run          = a cell array containg a vector of the framewise displacement for
%                         every timepoint for each run
%
%       snrmean         = mean temporal signal-to-noise ratio
%
%       snrmean_mid     = signal-to-noise ratio of the middle-slice of the timeseries
%
%       spatialsnr_mean = mean spatial signal-to-noise ratio
%
%       numBad          = a 1 x n vector containing the number of flagged timepoints for 
%                         each functional run (where n = number of functioanl runs)
%
% See also:

% Store output from previous steps into summary_info
summary_info        = cell(5, length(b.runs)+1);
summary_info(1,1)   = {b.curSubj};
summary_info(2:7,1) = {'Frame Displacement (mean)', 'Frame Displacement', 'Temporal SNR', 'Temporal SNR (mid)', 'Spatial SNR', 'Pct Bad Ts'};
for j = 1:length(b.runs)
    summary_info(2:7,j+1) = {FD_mean_run{j}, FD_run(j), snrmean(j), snrmean_mid(j), spatialsnr_mean(j), round(numBad(j))};
end

% Print summary info to screen
summary_info %#ok<NOPRT>

% Save summary information to QA directory
sum_name = fullfile(b.dataDir, b.QAdir, [b.curSubj '_QA_summary_information']);
save(sum_name, 'summary_info')

% Print all messages for this subject
b.messages = [b.messages 'END OF MESSAGES.\n'];
fprintf(b.messages);

% Clean up
clear summary_info
close all

end % summary_info
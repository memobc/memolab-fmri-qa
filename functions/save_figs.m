function [] = save_figs(fighandle, figname, b, figFormat, pdfReport, figwidth, figheight)
% SAVE_FIGS
% This function sets up basic plotting features for use with memolab_MRI_QA
%
% Note that if you're having problems with export_fig, you could switch
% this over to saveas, but it won't look as nice

% Use run length to specify figure height if figwidth and figheight are not
% given as inputs
if nargin < 6
    figwidth  = 6;
    figheight = figwidth .* length(b.runs) .* .25;
end

% Basic figure size & color properties
set(fighandle, 'Units', 'Inches', 'Position', [0, 0, figwidth, figheight])
set(fighandle, 'Color', 'w')

% Use export_fig to save figure
if pdfReport
    % Append to single QA_report pdf file
    export_fig(fighandle, fullfile(b.dataDir, b.QAdir, ['QA_report_',b.curSubj]), '-pdf', '-append');
    fprintf('Saving %s to file: %s\n', figname, fullfile(b.dataDir, b.QAdir, ['QA_report_' b.curSubj]));
else
    % Create a new image file for that figure
    export_fig(fighandle, fullfile(b.dataDir, b.QAdir, [figname '_' b.curSubj]), figFormat);
    fprintf('Saving %s to file: %s\n', figname, fullfile(b.dataDir, b.QAdir, [figname,'_',b.curSubj]));
end

% Clear figure handle
clear fighandle

end %save_figs

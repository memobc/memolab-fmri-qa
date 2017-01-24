function [b, numBad]= create_spike_regs(b, all_suspects, spikePrefix)
% CREATE_SPIKE_REGS  create spike regressors for bad timepoints
%
% This function uses frame displacement and global signal information to
% create spike regressors for bad timepoints. Spike regressors are combined
% with the original realignment parameters and written out to a file that
% can be pulled into any GLM (in place of the usual rp file).
%
%   USAGE:
%
%       [b, numBad] = create_spike_regs(b, all_suspects, spikePrefix)
%       
%       Input Arguments:
%
%       b            = memolab_MRI_batch structure
%
%       all_suspects = list of all bad timepoints as flagged by Art Repair
%
%       spikePrefix  = prefix for the spike regressors file. Note: can be
%                      blank, (i.e., '')
%
%       Output Arguments:
%
%       b       = memolab_MRI_batch structure
%
%       numBad = a double vector which contains the number of bad
%                timepoints flagged for each run
%
% See also: art_global_qa, memolab_batch_qa

for j = 1:length(b.runs)
    
    fprintf('\nWorking on %s\n',strcat(b.dataDir,b.runs{j}));
    
    % load original realignment parameters
    if size(b.rundir(j).rp,1) == 1
        motionFile = b.rundir(j).rp;
    else
        error('More than one motion file found');
    end
    motion = load(motionFile);
    
    % load ART suspects
    cursuspects    = all_suspects{j};
    bad_timepoints = unique(cursuspects);
    fprintf('-- Identified a total of %0.0f bad timepoints (%0.0f percent)\n', ...
        length(bad_timepoints),100*length(bad_timepoints)/size(motion,1));
    if (100*length(bad_timepoints)/size(motion,1)) > 20
        b.messages = [b.messages ...
            sprintf('Bad timepoints > 20 percent: consider excluding run %s.\n',b.runs{j})];
    end
    numBad(j) = (100*length(bad_timepoints)/size(motion,1));
    
    % write out list of bad timepoints
    filename = fullfile(b.dataDir, b.runs{j}, [spikePrefix 'bad_timepoints.txt']);
    fprintf('-- Writing out %s\n', filename);
    dlmwrite(filename, bad_timepoints, 'delimiter', '\t');
   
    % expand bad_timepoints into spike regressors
    if length(bad_timepoints) > 1
        spikereg = zeros(size(motion,1), length(bad_timepoints));
        for bad = 1:length(bad_timepoints)
            spikereg(bad_timepoints(bad),bad) = 1;
        end
    else
        spikereg = [];
    end
    
    % tack on spike regs to motion regressors
    allreg = [motion spikereg];
    
    % write out new rp regressors including spike regs
    filename = fullfile(b.dataDir, b.runs{j}, [spikePrefix 'spike_regs_rp.txt']);
    fprintf('-- Writing out %s\n\n', filename);
    dlmwrite(filename, allreg, 'delimiter', '\t');
    
    clear bad_timepoints spikereg rpparams allreg
end

end % create_spike_regs

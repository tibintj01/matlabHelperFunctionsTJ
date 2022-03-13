function [seqState, brainState] = assignTimesToBrainStates(stateInfo, timesToAssign)
% INPUTS:
% stateInfo is the struct returned from processHypnogram
% timesToAssign are times within the start time of the first state and the
% end time of the last state (e.g. theta cycle min times, spike times,
% etc.)
% OUTPUTS:
% seqState: the sequential state to which each time belongs (e.g. the first
% state in stateInfo or the 10th state, etc.)
% brainState: the brain state to which each time belongs (e.g. 0,1,2,3
% corresponding to Awk = 0, Unclass = 1, NREM = 2, REM = 3;)
%
% Omar Ahmed
% Version 1.0: 20180504

%%

% make sure you use interp1 (tested for 2017b, not earlier) and not earlier
% versions of interp1 or interp1_same_x_allowed (which is something 
seqState = interp1([stateInfo.startTime stateInfo(end).endTime], [1:(length([stateInfo.startTime])+1)], timesToAssign, 'previous'); %assign to the previous start time - works well and fast
brainState = NaN(size(seqState));
nonNan = find(isnan(seqState) == 0);
brainState(nonNan) = [stateInfo(seqState(nonNan)).state];
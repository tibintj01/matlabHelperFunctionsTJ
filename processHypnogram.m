function [allStateInfo] = processHypnogram(hyp, hypTimeAxis)
% Processes the hypnogram and returns a struct that correspond to each
% state transition
% Codes for Brain States (not needed here, just for information purposes)
% Awk = 0; Unclass = 1; NREM = 2;REM = 3;

% TO DO: Also return info about a "SLEEP EPOCH" which we can define as an epoch
% where there was no awake bout greater than X seconds, and X can be an
% input variable. By grouping and numbering these sleep epochs we can do
% sleep-epoch related analyses of spike rates and theta-gamma locking etc.

% NOTE:
% If hypTimeAxis = [0.5 1.0 1.5 ...], then the corresponding hyp values are
% for 0 - 0.5s, 0.5 - 1.0s, 1.0 - 1.5s, ...]
% Omar Ahmed, V1.0, 20180603

hypWindow = diff(hypTimeAxis(1:2));
hypOffset = hypWindow; % as the hypTimeAxis value corresponds to the end of the hypnogram calculation window

diffHyp = [0 diff(hyp)];
transitionInd = find(diffHyp ~= 0);
%transitionTimes = hypTimeAxis(transitionInd);

numTransitions = length(transitionInd);
numStates = numTransitions + 1;

allStateInfo = struct([]);
currStartInd = 1;
for i=1:numStates
    if i==numStates
        currEndInd = length(hyp);
    else
        currEndInd = transitionInd(i) - 1;
    end
    currState = unique(hyp(currStartInd:currEndInd));
    allStateInfo(i).state = currState;
    allStateInfo(i).startInd = currStartInd;
    allStateInfo(i).endInd = currEndInd;
    allStateInfo(i).startTime = hypTimeAxis(currStartInd) - hypOffset;
    allStateInfo(i).endTime = hypTimeAxis(currEndInd); % The hypTimeAxis aleady refers to the end time, so no offset
    allStateInfo(i).duration = allStateInfo(i).endTime - allStateInfo(i).startTime;
    allStateInfo(i).midTime = allStateInfo(i).startTime + (allStateInfo(i).duration/2); 
    if i==1
        allStateInfo(i).prevState = NaN;
    else
        allStateInfo(i).prevState = allStateInfo(i-1).state;
        allStateInfo(i-1).nextState = allStateInfo(i).state;
    end
    if i==numStates
        allStateInfo(i).nextState = NaN;
    else
        currStartInd = transitionInd(i); % set this for the next cycle
    end
end
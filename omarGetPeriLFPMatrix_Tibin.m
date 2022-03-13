function [eventLFP timeAxis err N trialsUsed] = omarGetPeriLFPMatrix_Tibin(lfp, Fs, eventBins, startTime, endTime, normAmp)
% helper function to average out the waveforms at the start,peak and min
% eventBins is in bins
% it is assumed that all the indices are the same length
% totalTime is in ms
% normTime, normAmp determine whether the time and amplitude are normalized
% from 0 to 1 or not
% - use INTERP1 when time averaging
% - amp averaging is easy
% XX - need to add error bars!!
% REMOVED normTIME for now

numEvents = length(eventBins);
relBins = [round(startTime*Fs) : round(endTime*Fs)];
numBins = length(relBins);
eventLFP = NaN(numEvents, numBins);
trialsUsed = [];
N = 0;
for i=1:numEvents
    currBins = eventBins(i) + relBins;
    if min(currBins) < 1 || max(currBins) > length(lfp)
        disp(['omarGetAveragedWaveforms: desired waveform exceeds LFP edges, trial ' num2str(i)]);
        continue;
    end
    if normAmp == 1
        currMin = min(lfp(currBins));
        currMax = max(lfp(currBins));
        currRange = currMax-currMin;
        scaledLFP = (lfp(currBins)-currMin)./currRange; % scales it from 0 to 1
        eventLFP(i,:) = scaledLFP;
    else
        eventLFP(i,:) = lfp(currBins);
    end
    trialsUsed = [trialsUsed; i];
    N = N + 1;
end
%avg = nanmean(eventLFP,1);
err = nanstd(eventLFP, [], 1) ./ sqrt(N);
timeAxis = relBins ./ Fs;

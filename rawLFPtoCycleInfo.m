function [cycleInfo,extractedTimeAxis,extractedLFP] = rawLFPtoCycleInfo(timeAxis,lfpValues,lowFreq,highFreq)
%
    Fs=1/median(diff(timeAxis));
    %filterOrder=2;
    %zscore=0;
    %[filteredLFP]=filterLFP(lfpValues,Fs,lowFreq,highFreq,filterOrder,zscore);
    %saveDir='cycleProcessedData'

    %minTime=min(timeAxis);
    %maxTime=max(timeAxis);

    [extractedTimeAxis, extractedLFP, cycleInfo] = ...
        getCycleInfoLFP(timeAxis, filteredLFP, Fs, lowFreq, highFreq, 1);


end


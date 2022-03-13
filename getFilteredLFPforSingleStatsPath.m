function [decLFP]=getFilteredLFPforSingleStatsPath(currChSingleCycleFilePath,lowFastFreq,highFastFreq)

 disp(sprintf('extracting cycles from ch %d........',i))
        %load single cycle info
        data=load(currChSingleCycleFilePath);

        %extract raw decimated LFP corresponding to current channel
        fileName=getFileNameFromPath(currChSingleCycleFilePath);
        ch=str2num(fileName((end-5):(end-4)));
        [decLFP,decFs]=getDecLFPforCh(ch);

                %get filtered LFP in desired (broad) frequency band
                %[lowFastFreq,highFastFreq]=getFastFreqBand(freqBandName);
                %decLFP=filterLFP(decLFP,decFs,data.lowFreq,data.highFreq,2,1);
                decLFP=filterLFP(decLFP,decFs,lowFastFreq,highFastFreq,2,1);

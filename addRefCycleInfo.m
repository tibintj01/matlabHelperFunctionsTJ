close all; clear all; clc
dataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles';

filePaths=getFilePathsRegex(dataDir,'*mat');

for fi=1:length(filePaths)
    %for fi=75:length(filePaths)
    fi
    currFilePath=filePaths{fi};
    currFileName=getFileNameFromPath(currFilePath);
    fileBaseName=currFileName(1:(end-4));
    
    [currSessionNum,currSessName]=getSessionNum(fileBaseName);

    data=load(currFilePath);
    
    [refThetaCh] = getRefThetaChNumForSession(fileBaseName);
    if(isnan(refThetaCh))
        continue
    end
    refThetaChCyclePath=sprintf('./hc11ThetaCycleData/%s_Ch%d_6-12_minmax.mat',currSessName,refThetaCh);
    
    refCycleInfo=load(refThetaChCyclePath);
    
    refMinTimes=refCycleInfo.cycleMinTimes;
    refMaxTimes=refCycleInfo.cycleMaxTimes;
    
    for di=1:2
        if(di==1)
            currDirStr='right';
            currSpikeTimes=data.rightSpikeTimes;
        else
            currDirStr='left';
            currSpikeTimes=data.leftSpikeTimes;
        end
    
             [spikeAssignMax spikeAssignMin spikeAssignZero] = assignSpikesToCycles2017(currSpikeTimes, refMinTimes, refMaxTimes);

                mode=2; %phase from max to max
                %mode=1;%phase from startMin to endMin
                calcTimeOffsets=1;
                [refThetaSpikePhases spikeOffsets spikeOffsetsEnd] = assignPhaseToSpikes(currSpikeTimes, spikeAssignMax, spikeAssignMin, refMinTimes, refMaxTimes, mode, calcTimeOffsets, spikeAssignZero);
            
         spikeCycleIDs=spikeAssignMin;
         isFirstSpikeInCycle=[diff(spikeCycleIDs(:)); NaN]>0;
         refThetaSpikePhasesOnlyCycleFirsts=refThetaSpikePhases(isFirstSpikeInCycle);
         refThetaSpikeTimesOnlyCycleFirsts=currSpikeTimes(isFirstSpikeInCycle);
                
        if(di==1)
            rightSpikesRefPhasesMaxToMax=refThetaSpikePhases;
            rightSpikeRefCycleIDs=spikeCycleIDs;
            rightIsFirstSpikeInRefCycle=isFirstSpikeInCycle;
            rightSpikePhasesOnlyRefCycleFirsts=refThetaSpikePhasesOnlyCycleFirsts;
            rightSpikeTimesOnlyRefCycleFirsts=refThetaSpikeTimesOnlyCycleFirsts;
            save(currFilePath,'refThetaCh', 'refThetaChCyclePath','rightSpikesRefPhasesMaxToMax','rightSpikeRefCycleIDs','rightIsFirstSpikeInRefCycle','rightSpikePhasesOnlyRefCycleFirsts','rightSpikeTimesOnlyRefCycleFirsts','-append')
        else
            leftSpikesRefPhasesMaxToMax=refThetaSpikePhases;
            leftSpikeRefCycleIDs=spikeCycleIDs;
            leftIsFirstSpikeInRefCycle=isFirstSpikeInCycle;
            leftSpikePhasesOnlyRefCycleFirsts=refThetaSpikePhasesOnlyCycleFirsts;
            leftSpikeTimesOnlyRefCycleFirsts=refThetaSpikeTimesOnlyCycleFirsts;
            save(currFilePath,'refThetaCh', 'refThetaChCyclePath','leftSpikesRefPhasesMaxToMax','leftSpikeRefCycleIDs','leftIsFirstSpikeInRefCycle','leftSpikePhasesOnlyRefCycleFirsts','leftSpikeTimesOnlyRefCycleFirsts','-append')
        end
    end
    
    
       
end
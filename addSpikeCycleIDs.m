
close all; clear all; clc
%dataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles';

processedDataDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData';
unitDataDir=fullfile(processedDataDir,'unitSpikeInfoStructs');
dataDir=unitDataDir;


filePaths=getFilePathsRegex(dataDir,'*mat');


startFileID=1;

for fi=startFileID:length(filePaths)
       fi
    currFilePath=filePaths{fi};
    
    currFileName=getFileNameFromPath(currFilePath);
    fileBaseName=currFileName(1:(end-4));
    
   
    
    data=load(currFilePath);
    
    cycleData=load(data.refChCycleInfoPath);
    cycleStartTimes=cycleData.cycleMaxTimes;
    
    posInfo=load(data.unitInfo.positionInfoForSessionSavePath);
    
    spikeTimeDirectionAssignment=data.unitInfo.spikeTimeDirectionAssignment;
    
    leftSpikeCycleIDs=NaN;
    rightSpikeCycleIDs=NaN;
    
    for di=1:2
        unitDirSpikeTimes=data.unitInfo.spikeTimes(spikeTimeDirectionAssignment==di);
        
        if(di==2)
            leftSpikeCycleIDs=interp1(cycleStartTimes,1:length(cycleStartTimes),unitDirSpikeTimes);
        else
            rightSpikeCycleIDs=interp1(cycleStartTimes,1:length(cycleStartTimes),unitDirSpikeTimes);
        end
        
        
    end
    
    save(currFilePath,'leftSpikeCycleIDs', 'rightSpikeCycleIDs','-append')
   
    
end
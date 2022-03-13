%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load unit struct file names
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; clc
processedDataDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData';
unitDataDir=fullfile(processedDataDir,'unitSpikeInfoStructs');


filePaths=getFilePathsRegex(unitDataDir,'*mat');

startFile=1;
for fi=startFile:length(filePaths)
    %for fi=76:length(filePaths)
       fi
    currFilePath=filePaths{fi};
    currFileName=getFileNameFromPath(currFilePath);
    fileBaseName=currFileName(1:(end-4));
    
    dataStruct=load(currFilePath);
    data=dataStruct.unitInfo;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %load gamma cycle timing data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    maxRippleChThisShank=data.maxRippleChThisShank;
    lowFreq=40;
    highFreq=80;
    cycleInfoFileName=fullfile(saveLFPcycleInfoDir,sprintf('%s_Ch%d_%d-%d_minmax.mat',currSessionName,maxRippleChThisShank,lowFreq,highFreq));
    currUnitCycleInfo=load(cycleInfoFileName);

    minTimes=currUnitCycleInfo.cycleMinTimes;
    maxTimes=currUnitCycleInfo.cycleMaxTimes;


    numCycles=length(currUnitCycleInfo.cycleMaxTimes);
    [spikeAssignMax spikeAssignMin spikeAssignZero] = assignSpikesToCycles2017(unitSpikeTimes, minTimes, maxTimes);

    mode=2; %phase from max to max
    %mode=1;%phase from startMin to endMin
    calcTimeOffsets=1;
    [spikePhase spikeOffsets spikeOffsetsEnd] = assignPhaseToSpikes(unitSpikeTimes, spikeAssignMax, spikeAssignMin, minTimes, maxTimes, mode, calcTimeOffsets, spikeAssignZero);

    posInfo=load(data.positionInfoForSessionSavePath);
    unitSpikeGammaPhases=spikePhase(:);
    
    
  speedPerTimeStep=posInfo.speedAlongTrackPerTimeMSec;
  speedTimeAxis=posInfo.positionTimeAxisSec;
    
    for di=1:2
        
      
        
        if(di==2)
            currFieldDirection='leftward'; 
            goingBackwardPerTimeStep=speedPerTimeStep>0; 
        elseif(di==1)
            currFieldDirection='rightward'; 
            goingBackwardPerTimeStep=speedPerTimeStep<0;
        end
        
         currDirSpikeTimes=data.spikeTimes(data.spikeTimeDirectionAssignment==di);
         
         currDirSpikeGammaPhases=unitSpikeGammaPhases(data.spikeTimeDirectionAssignment==di);
         currDirSpikeRunningSpeeds=interp1(speedTimeAxis,speedPerTimeStep,currDirSpikeTimes);
        
        
        
        
    end

end

close all; clear all; clc
dataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles';

%GAUSSKERNELWIDTHSEC=0.1;
GAUSSKERNELWIDTHSEC=0.01;

filePaths=getFilePathsRegex(dataDir,'*mat');

for fi=1:length(filePaths)
       fi
    currFilePath=filePaths{fi};
    currFileName=getFileNameFromPath(currFilePath);
    fileBaseName=currFileName(1:(end-4));
    
    currSessionNum=getSessionNum(fileBaseName);

    data=load(currFilePath);
     
     %currUnitInfo=data.unitInfo;
    figure; histogram(data.unitInfo.spikePhases)
     %localCycleInfo=load(currUnitInfo.cycleInfoFileName);
     

    cycleData=load(data.unitInfo.cycleInfoFileName);
    %cycleStartTimes=cycleData.cycleMaxTimes(1:(end-1));
    %cycleEndTimes=cycleData.cycleMaxTimes(2:end);
    cycleStartTimes=cycleData.cycleMinTimes(1:(end-1));
    cycleEndTimes=cycleData.cycleMinTimes(2:end);
    numCycles=length(cycleStartTimes);
    figure; 
    for di=1:2
        if(di==1)
            dirStr='right';
        else
            dirStr='left';
        end
        
        currSpikeTimes=data.(sprintf('%sSpikeTimes',dirStr));
        %currSpikePhases=data.(sprintf('%sSpikePhases',dirStr));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %PHASE WRT ONE WELL CHOSEN REFERENCE CHANNEL
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        currSpikePhases=data.(sprintf('%sSpikePhases',dirStr));
        %{
        subplot(2,1,di)
        histogram(currSpikePhases,100)
        %}
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %loop through each cycle
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for ci=1:numCycles
            currCycleStartTime=cycleStartTimes(ci);
            currCycleEndTime=cycleEndTimes(ci);
            
            spikesInCurrCycleIdxes=currSpikeTimes>=currCycleStartTime & ...
                currSpikeTimes<currCycleEndTime;
            
            currCycleSpikePhases=currSpikePhases(spikesInCurrCycleIdxes);
            currCycleSpikeTimes=currSpikeTimes(spikesInCurrCycleIdxes);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %identify in cycle spike rank (min to min cycle) of each spike time
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %{
            if(length(currCycleSpikeTimes)>3)
                figure
                plot(currCycleSpikeTimes,currCycleSpikePhases,'k.','MarkerSize',30)
                hold on
                ylim([0 360])
                plot([currCycleStartTime currCycleStartTime],ylim,'k--','LineWidth',4)
                plot([currCycleEndTime currCycleEndTime],ylim,'k--','LineWidth',4)

                  disp('')
                  close all
            end
            %}
            
        end
        
       
    end
     close all
     
    %unitCycleCouplingInfo.
    %save(currFilePath,'rightSpikeInCycleRanks','leftSpikeInCycleRanks','-append') 
end
%}


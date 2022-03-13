close all; clear all; clc
dataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles';
cycleDataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/';
showPlots=1;
zScorePlot=1;
%zScorePlot=0;

minSpeed=0.05;

cycleHalf=0.06; %60 msec
cycleDur=0.12;


maxSpaceDiff=1;
maxSpaceDiff=0.5;
%maxSpaceDiff=0.3;
maxSpaceDiff=1;
spaceFracMax=0.3*0.3;
spaceFracMax=0.15;
spaceFracMax=0.5;
spaceFracMax=0.2;
%spaceFracMax=0.15;
%spaceFracMin=0.05;
spaceFracMin=0;
spaceFracMin=0.01;
spaceFracMin=0;
%spaceFracMin=0.05;
%spaceFracMax=0.05;

%{
minSpeedForThetaDiff=0.2;
maxSpeedForThetaDiff=0.6;
%}


minSpeedForThetaDiff=0;
maxSpeedForThetaDiff=10;


lastFrac=0.35;
lastFrac=0.33333;
%lastFrac=0.4;
lastFrac=0.5;
%lastFrac=0.4;

%lastFrac=0.2
cofiringThresh=0.2;
cofiringThresh=0.1;
%cofiringThresh=0.08;

lateCycleTimeThresh=cycleDur*(1-lastFrac);
earlyCycleTimeThresh=cycleDur*lastFrac;
numBins=30;
numBins=50;
%{
lateCycleTimeThresh=cycleDur*0.6;
earlyCycleTimeThresh=cycleDur*0.4;
%}

%{
lateCycleTimeThresh=cycleDur*0.5;
earlyCycleTimeThresh=cycleDur*0.5;
%}

%Cicero_09172014 has only 1 strongly precessing cell used
%sessionNames={'Achilles_10252013','Achilles_11012013','Buddy_06272013', 'Cicero_09102014','Cicero_09172014','Gatsby_08282013'};
sessionNames={'Achilles_10252013','Achilles_11012013','Buddy_06272013', 'Cicero_09102014','Gatsby_08282013'};

%sessionNames={'Achilles_11012013'};
%sessionNames={'Gatsby_08282013'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%loop through each session and load per cycle spike time data across units
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for si=1:length(sessionNames) 
    currSessName=sessionNames{si};
    currSeqDataFilePath=sprintf('perUnitPerCycleTimingData%s.mat',currSessName);
    
    currThetaSeqData=load(currSeqDataFilePath);
     numUnits=size(currThetaSeqData.avgSpikeTimesByCycleAndUnit,2);
     numCycles=currThetaSeqData.numCycles;
     sortedMeanCenterIdxPerDi=currThetaSeqData.sortedMeanCenterIdxPerDi;
     timeOfCycleTrough=currThetaSeqData.timeOfCycleTrough;
     
     numCyclesCofiring=zeros(numUnits,numUnits,2);
       
     [refUnitDataLeftward,refUnitDataRightward] = getRefUnitData(currSessName)
     

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %for each direction, loop through each pair of units and each cycle 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for di=1:2
        
        if(di==1)
            currDiStr='rightward';
            
            if(strContainsCircSessionName(currSessName))
                continue
            end
            refUnitData=refUnitDataRightward;
        else
            currDiStr='leftward';
            refUnitData=refUnitDataLeftward;
        end
        
        timeDifferencesPerCycle=NaN(numUnits,numUnits,numCycles);
        currTimeInMetaFieldPerCycle=NaN(numUnits,numUnits,numCycles);
        numCyclesThreshPerRow=NaN(numUnits,1);
        numCyclesMaxPerRow=NaN(numUnits,1);
        alreadyProcessed=false(numUnits,numUnits);
        
        currFieldCentersPerCycle=currThetaSeqData.unitCenterPerCycleAndField(:,:,di);
        currFieldStartTimesPerCycle=currThetaSeqData.unitEntryTimePerCycleAndField(:,:,di);
        currFieldEndTimesPerCycle=currThetaSeqData.unitExitTimePerCycleAndField(:,:,di);
        
        
        avgSpikeTimesByCycleAndUnit=squeeze(currThetaSeqData.avgSpikeTimesByCycleAndUnit(:,:,di));
        
        numCycles=size(avgSpikeTimesByCycleAndUnit,1);
        numCells=size(avgSpikeTimesByCycleAndUnit,2);
        
        minPosition=min(refUnitData.positionPerTimeStep);
        maxPosition=max(refUnitData.positionPerTimeStep);
        totalLength=maxPosition-minPosition;
        
        %loop through cycles, find all active cells in current cycle
        numSeqs=0;
        maxRank=7;
        positionColors=jet(100);
        maxNumCyclesDisp=50;
         maxNumCyclesDisp=10;
        for ci=1:numCycles
            currCycleAllUnitAvgSpikeTimes=avgSpikeTimesByCycleAndUnit(ci,:);
            currCycleAllUnitCentersM=currFieldCentersPerCycle(ci,:);
            
            currCycleNonNanSpikeTimes=currCycleAllUnitAvgSpikeTimes(~isnan(currCycleAllUnitAvgSpikeTimes));
            currCycleNonNanAllUnitCentersM=currCycleAllUnitCentersM(~isnan(currCycleAllUnitAvgSpikeTimes));
            
            cycleSpan=range(currCycleNonNanSpikeTimes);
            %{
            if(length(currCycleNonNanSpikeTimes)<5 || cycleSpan<0.06)
                continue
                end
            %}
            
            
            if(length(currCycleNonNanSpikeTimes)<3)
                continue
            end
            numSeqs=numSeqs+1;
            if(numSeqs>maxNumCyclesDisp)
                numSeqs=1;
            end
           
            [sortedSpikeTimes,sortIdx]=sort(currCycleNonNanSpikeTimes);
            currCycleNonNanAllUnitCentersMsorted=currCycleNonNanAllUnitCentersM(sortIdx);
            
            minPositionInSeq=min(currCycleNonNanAllUnitCentersMsorted);
            currCycleNonNanAllUnitCentersMsorted=currCycleNonNanAllUnitCentersMsorted-minPositionInSeq;
            
            totalPosSpanOfSeq=range(currCycleNonNanAllUnitCentersMsorted);
            for si=1:length(sortedSpikeTimes)
                %rankID=min(si,maxRank);
                
                %position
                
                %plotRasterStyle(currCycleNonNanSpikeTimes(si),numSeqs,NaN,NaN,rankColors(rankID,:))
                positionColorID=ceil((currCycleNonNanAllUnitCentersMsorted(si)+0.001)/totalPosSpanOfSeq*100);
                
                positionColorID=min(100,positionColorID);
                %plotRasterStyle(currCycleNonNanSpikeTimes(si),numSeqs,NaN,NaN,positionColors(positionColorID,:),1)
                plotRasterStyle(currCycleNonNanSpikeTimes(si),numSeqs,NaN,NaN,NaN,1)

                hold on
            end
            xlim([0 .1])
            ylim([0 maxNumCyclesDisp])
            %colormap(jet)
            %cb=colorbar;
            %ylabel(cb,'rank in sequence')
            drawnow
        end
        
    end
    
end


close all; clear all; clc
%dataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles';
%HC3 data
processedDataDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData';
unitDataDir=fullfile(processedDataDir,'unitSpikeInfoStructs');
dataDir=unitDataDir;

getLocalGaussSpikeRate=0;
getLocalGaussSpikeRate=1; %needed for estimating phase per cycle and time point


spaceTimePlotSaveDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3Images/worldLineVsPhase';
touchDir(spaceTimePlotSaveDir)

showPlots=0;

filePaths=getFilePathsRegex(dataDir,'*mat');
totalNumFields=0;

resetGaussSpikeWidth=1;

resetGaussSpikeWidth=0;
  %kernelWidthSec=0.1; %sec
    kernelWidthSec=0.01; %sec
    
useManualTimeFieldBoundForNorm=1;

startFileIdx=1;
%startFileIdx=672;
%startFileIdx=650;
for i=startFileIdx:length(filePaths)
%for i=1:1
    
    %try
        
    i/length(filePaths)
    currDataFilePath=filePaths{i};
    
    currFileName=getFileNameFromPath(currDataFilePath);
    fileBaseName=currFileName(1:(end-4));
    
    %{
    if(~strcmp(debugTimeFilePath,currDataFilePath))
       %continue
    end
    %}
    
    currData=load(currDataFilePath);
    
    if(isfield(currData,'directionSpecificStats'))
        if(~isempty(currData.directionSpecificStats))
         %continue
        end
    end
    
   
    %{
    if(max(currData.numFieldsPerDir)<2)
        continue
    end
    if(isfield(currData,'directionSpecificStats'))
        %continue
    end
    %}
    
    if(~isfield(currData,'entryTimePerDirPerFieldPerLap'))
        continue
    end
    
    entryTimePerDirPerFieldPerLap=currData.entryTimePerDirPerFieldPerLap;
    exitTimePerDirPerFieldPerLap=currData.exitTimePerDirPerFieldPerLap;
    
    entryPosPerDirPerFieldPerLap=currData.entryPosPerDirPerFieldPerLap;
    exitPosPerDirPerFieldPerLap=currData.exitPosPerDirPerFieldPerLap;
    
    entryTimeIdxesPerDirPerFieldPerLap=currData.entryTimeIdxesPerDirPerFieldPerLap;
    exitTimeIdxesPerDirPerFieldPerLap=currData.exitTimeIdxesPerDirPerFieldPerLap;
    
    availableDirections=1:2;
    
    if(length(currData.leftwardFieldStartEndM)==0)
        availableDirections=1;
    end
    if(length(currData.rightwardFieldStartEndM)==0)
        availableDirections=2;
    end
  
      directionSpecificStats=[];
    for di=availableDirections
        if(showPlots)
            close all
              figure
          
        end
        currFieldDirection='';
        spikeTimeDirectionAssignment=currData.unitInfo.spikeTimeDirectionAssignment;
        
        if(di==2)
            currFieldDirection='leftward';  
            if(showPlots)
                
                imshow(currData.unitInfo.saveImgPathLeftward)

            end
            
            directionalSpikeTimes=currData.unitInfo.spikeTimes(spikeTimeDirectionAssignment==2);
            directionalSpikePhases=currData.unitInfo.leftSpikePhases;
            directionalSpikePositions=currData.unitInfo.leftSpikePositions;
            
        elseif(di==1)
            currFieldDirection='rightward';  
            if(showPlots)
               imshow(currData.unitInfo.saveImgPathRightward)
            end
            directionalSpikeTimes=currData.unitInfo.spikeTimes(spikeTimeDirectionAssignment==1);
            directionalSpikePhases=currData.unitInfo.rightSpikePhases;
            directionalSpikePositions=currData.unitInfo.rightSpikePositions;
        end
        
        if(showPlots)
         autoArrangeFigures(1,2,1);
        end
        if(~(isfield(currData,sprintf('%sFieldStartEndM',currFieldDirection))))
            continue
        end
        
        fieldStartEndPosFieldName=sprintf('%sFieldStartEndM',currFieldDirection);

       startEndPositions=currData.(fieldStartEndPosFieldName);
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %in HC3 manual field start and end bounds are all in one array
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       numBounds=length(startEndPositions);
       %numFields=numBounds/2;
       numFields=size(entryTimePerDirPerFieldPerLap,2);
       
       
       startPositions=NaN(numFields,1);
       endPositions=NaN(numFields,1);
       fieldCount=1;
       for bi=1:numBounds
           if(mod(bi,2)==1)
            startPositions(fieldCount)=startEndPositions(bi);
           else
            endPositions(fieldCount)=startEndPositions(bi);
            fieldCount=fieldCount+1;
           end
       end
   
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %cycle data from pyr layer theta
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       cycleData=load(currData.refChCycleInfoPath);
       cycleStartTimes=cycleData.cycleMaxTimes;
       
        posInfo=load(currData.unitInfo.positionInfoForSessionSavePath);
        positionTimeAxis=posInfo.positionTimeAxisSec;
       
       
       cycleStartTimes=cycleStartTimes(cycleStartTimes>=min(positionTimeAxis) & cycleStartTimes<=max(positionTimeAxis));
       
       allSpikeTimesAndPhases=[directionalSpikeTimes(:) directionalSpikePhases(:)];
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %local spike phase based on gaussian kernel peak per theta cycle
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %cycleStartTimes=
     
      
       %[localCycleMeanSpikePhasePerTime, localCycleMeanSpikeTimePerTime]=getLocalCycleMeanPhasePerTime(allSpikeTimesAndPhases,positionTimeAxis,cycleStartTimes);
     
   if(getLocalGaussSpikeRate)
       %if(~resetGaussSpikeWidth && isfield(currData,directionSpecificStats) && isfield(currData.directionSpecificStats,currFieldDirection) && isfield(currData.directionSpecificStats.(currFieldDirection),'localSpikePhasePerTime'))
       if(~resetGaussSpikeWidth && isfield(currData,'directionSpecificStats') && isfield(currData.directionSpecificStats,currFieldDirection) && isfield(currData.directionSpecificStats.(currFieldDirection),'localSpikePhasePerTime'))

           localSpikePhasePerTime=currData.directionSpecificStats.(currFieldDirection).localSpikePhasePerTime;
           gaussTimeAxis=currData.directionSpecificStats.(currFieldDirection).gaussTimeAxis;
           localGaussSpikeRatePerTime=currData.directionSpecificStats.(currFieldDirection).localGaussSpikeRatePerTime;
           if(isfield(currData,'rightLocalGaussSpikeRatePerTime'))
            rightLocalGaussSpikeRatePerTime=currData.rightLocalGaussSpikeRatePerTime;
           else
               rightLocalGaussSpikeRatePerTime=NaN(size(gaussTimeAxis));
           end
            if(isfield(currData,'leftLocalGaussSpikeRatePerTime'))
             leftLocalGaussSpikeRatePerTime=currData.leftLocalGaussSpikeRatePerTime;
           else
               leftLocalGaussSpikeRatePerTime=NaN(size(gaussTimeAxis));
            end
          
            if(isfield(currData,'leftLocalGaussSpikeRatePerTime'))
                gaussTimeAxisLR=currData.gaussTimeAxis;
            else
                gaussTimeAxisLR=NaN(size(gaussTimeAxis));
            end
           
            decFactor=1;
       else
          [gaussTimeAxis,localGaussSpikeRatePerTime]=localGaussRateSpikeTrain(directionalSpikeTimes(:)-min(positionTimeAxis), kernelWidthSec,0,max(positionTimeAxis)-min(positionTimeAxis));
             gaussTimeAxis=gaussTimeAxis+min(positionTimeAxis);


           localSpikePhasePerTime=getPeakPhasePerThetaCyclePerTime(gaussTimeAxis,localGaussSpikeRatePerTime,cycleStartTimes,positionTimeAxis);
             decFactor=2;
       end
       
    end
       
       numFields=size(entryTimePerDirPerFieldPerLap,2);
       numLaps=size(entryTimePerDirPerFieldPerLap,3);
       
       fieldEntryTimePerTime=NaN(size(positionTimeAxis));
       fieldExitTimePerTime=NaN(size(positionTimeAxis));
       fieldNumPerTime=NaN(size(positionTimeAxis));
       lapNumberPerTime=NaN(size(positionTimeAxis));
       
       fieldEntryPosPerTime=NaN(size(positionTimeAxis));
       fieldExitPosPerTime=NaN(size(positionTimeAxis));
       
       avgFieldDurations=NaN(numFields,1);
       avgFieldWidths=NaN(numFields,1);
       
       fieldWidthsPerLap=NaN(numFields,numLaps);
       fieldDurationsPerLap=NaN(numFields,numLaps);
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %put down field characteristics into time domain
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       for fi=1:numFields
               currFieldDurations=squeeze(exitTimePerDirPerFieldPerLap(di,fi,:)-entryTimePerDirPerFieldPerLap(di,fi,:));

           %avgFieldDurations(fi)=nanmean(deleteoutliers(currFieldDurations));
           fieldDurationsPerLap(fi,:)=currFieldDurations;
           
           currFieldWidths=squeeze(exitPosPerDirPerFieldPerLap(di,fi,:)-entryPosPerDirPerFieldPerLap(di,fi,:));
           %avgFieldWidths(fi)=nanmean(deleteoutliers(currFieldWidths));
           fieldWidthsPerLap(fi,:)=abs(currFieldWidths);
           
           for li=1:numLaps
             currEntryIdx=entryTimeIdxesPerDirPerFieldPerLap(di,fi,li);
             currExitIdx=exitTimeIdxesPerDirPerFieldPerLap(di,fi,li);
             
            
             if(currEntryIdx>0)
                 fieldNumPerTime(currEntryIdx:currExitIdx)=fi;
                  lapNumberPerTime(currEntryIdx:currExitIdx)=li;
                  
                 fieldEntryTimePerTime(currEntryIdx:currExitIdx)=entryTimePerDirPerFieldPerLap(di,fi,li);
                 fieldExitTimePerTime(currEntryIdx:currExitIdx)=exitTimePerDirPerFieldPerLap(di,fi,li);

                 fieldEntryPosPerTime(currEntryIdx:currExitIdx)=entryPosPerDirPerFieldPerLap(di,fi,li);
                 fieldExitPosPerTime(currEntryIdx:currExitIdx)=exitPosPerDirPerFieldPerLap(di,fi,li);
             end
             
           end  
       end
       
       if(max(abs(fieldWidthsPerLap(:)))==0)
           continue
       end
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %create flag array indicating laps where field length is within 15 percent of expected and
       %field duration is greater than expected length times max speed 
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %{
       if(di==1)
           %expectedWidthPerField=currData.(fieldEndPosFieldName)-currData.(fieldStartPosFieldName);
           expectedWidthPerField=endPositions-startPositions;
       elseif(di==2)
           expectedWidthPerField=-(currData.(fieldStartPosFieldName)-currData.(fieldEndPosFieldName));
       end
       %}
       
        expectedWidthPerField=abs(endPositions-startPositions);
       
       %{
       if(expectedWidthPerField<0)
           disp('') 
     end
       %}
       
       
       expectedWidthPerField=expectedWidthPerField(:);
       expectedWidthPerField=repmat(expectedWidthPerField,1,size(fieldWidthsPerLap,2));
       
       %{
       removeWidthRows=[];
       for wr=1:size(expectedWidthPerField,1)
           if(sum(isnan(expectedWidthPerField(wr,:)))==size(expectedWidthPerField,2))
               removeWidthRows=[removeWidthRows wr];
           end
       end
       expectedWidthPerField(removeWidthRows,:)=[];
       %}
       
       
       %thresholdWidthDeviations=expectedWidthPerField*0.15;
        thresholdWidthDeviations=0.05; %within 5 cm of expected width
        thresholdWidthDeviations=0.1; %within 10 cm of expected width
       
         %maxSpeed=2;
       maxSpeed=1.3;
       minFieldDurations=abs(expectedWidthPerField)/maxSpeed;
       
       goodFieldWidthLaps=(abs(fieldWidthsPerLap-expectedWidthPerField)<=thresholdWidthDeviations);
       goodFieldDurationLaps=fieldDurationsPerLap>minFieldDurations;
       
        isGoodLapForEachField=goodFieldWidthLaps & goodFieldDurationLaps;
        
        if(max(isGoodLapForEachField(:))==0)
            continue
        end
        
         for fi=1:numFields
             goodLapsForField=isGoodLapForEachField(fi,:);

             goodFieldWidths=fieldWidthsPerLap(fi,goodLapsForField);
             if(~isempty(goodFieldWidths))
                avgFieldWidths(fi)=nanmean(deleteoutliers(goodFieldWidths));
                 %avgFieldWidths(fi)=nanmean((goodFieldWidths));
             else
                avgFieldWidths(fi)=NaN;
             end
             
             goodFieldDurs=fieldDurationsPerLap(fi,goodLapsForField);
             %goodFieldDurs(goodFieldDurs>(abs(avgFieldWidths(fi))/0.1))=NaN;%field width divided by min speed of 0.1m/sec
             [durN,durEdges,durBins]=histcounts(goodFieldDurs,30);
             [~,mostCommonDurBin]=max(durN);
             
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             %TIME NORMALIZATION
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             %typicalDur=nanmean(goodFieldDurs(durBins==mostCommonDurBin))*0.975;
             typicalDur=nanmean(goodFieldDurs(durBins==mostCommonDurBin));
             maxDurForAvg=typicalDur*2;
             %goodFieldDurs(goodFieldDurs>(abs(avgFieldWidths(fi))/0.05))=NaN;%field width divided by min speed of 0.1m/sec
          goodFieldDurs(goodFieldDurs>maxDurForAvg)=NaN;
             %abs(avgFieldWidths(fi)/0.1)
             if(~isempty(goodFieldDurs))
                avgFieldDurations(fi)=nanmean(deleteoutliers(goodFieldDurs));
                typicalDistance=nanmean(abs(typicalDur-goodFieldDurs));
                if(~isnan(typicalDur) && typicalDistance<0.25)
                    avgFieldDurations(fi)=typicalDur;
                end
                
               
             else
                 avgFieldDurations(fi)=NaN;
             end
             %figure; histogram(goodFieldDurs,30)
         end
        
         avgFieldWidths
         avgFieldDurations
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %compute time and distance within local field for each time point
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %localDistInFieldPerTime=currData.positionPerTimeStep;
        %localTimeInFieldPerTime=currData.positionTimeAxis;
        positionPerTimeStep=posInfo.posAlongTrackPerTimeM;
        positionTimeAxis=posInfo.positionTimeAxisSec;
        
        localDistInFieldPerTime=NaN(size(positionPerTimeStep));
        localTimeInFieldPerTime=NaN(size(positionTimeAxis));
        
        localNormDistInFieldPerTime=NaN(size(localDistInFieldPerTime));
        localNormTimeInFieldPerTime=NaN(size(localTimeInFieldPerTime));
        
        numTimePts=length(localDistInFieldPerTime);
        
        isCircMaze=currData.isCircMaze;
        for ti=1:numTimePts
            
            if(~isnan(fieldEntryPosPerTime(ti)) && isnan(fieldNumPerTime(ti)))
                disp('')
                
            end
            currFieldNum=fieldNumPerTime(ti);
            if(isnan(currFieldNum))
                %localDistInFieldPerTime(ti)=NaN;
                %localTimeInFieldPerTime(ti)=NaN;
                continue
            end

            currLapNum=lapNumberPerTime(ti);
            if(~isGoodLapForEachField(currFieldNum,currLapNum))
                %localDistInFieldPerTime(ti)=NaN;
                %localTimeInFieldPerTime(ti)=NaN;
                continue
            end
            

          
             localDistInFieldPerTime(ti)=abs(positionPerTimeStep(ti)-fieldEntryPosPerTime(ti));
             
               %number of cycles elapsed in field (count of troughs)
            %cycleStartTimes>=fieldEntryTimePerTime(ti))
            [~,currCycleNum]=min(abs(cycleStartTimes-positionTimeAxis(ti)));
            
            if(cycleStartTimes(currCycleNum)>positionTimeAxis(ti))
                currCycleNum=currCycleNum-1;
            end
            
            [~,entryCycleNum]=min(abs(cycleStartTimes-fieldEntryTimePerTime(ti)));
            if(cycleStartTimes(entryCycleNum)>fieldEntryTimePerTime(ti))
                entryCycleNum=entryCycleNum-1;
            end
            
             %localTimeInFieldPerTime(ti)=abs(positionTimeAxis(ti)-fieldEntryTimePerTime(ti));
             localTimeInFieldPerTime(ti)=currCycleNum-entryCycleNum;
            
            currFieldAvgWidth=avgFieldWidths(currFieldNum);
            localNormDistInFieldPerTime(ti)=localDistInFieldPerTime(ti)/abs(currFieldAvgWidth);
            
            if(useManualTimeFieldBoundForNorm)
                %currFieldNum reset for each time point
                %thetaBasedTimeInfo=currData.thetaBasedTimeVsPhaseInfo.(currFieldDirection).sprintf('field%d',currFieldNum);
                timeBoundTable=currData.cycleTimeFieldBoundPerDirPerField;
                
                %currFieldManualDuration=currData.(sprintf('manualTimeField%dEnd%sSec',currFieldNum,currFieldDirection));
                currFieldManualDuration=timeBoundTable.(currFieldDirection)(fi);
                
                %currFieldAvgDur=currFieldManualDuration;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %MOST RELEVANT TIME SCALE?
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 %currFieldAvgDur=currFieldManualDuration*0.8; 
                 currFieldAvgDur=currFieldManualDuration;  
                 %currFieldAvgDur=currFieldManualDuration*0.75; 
                %sprintf('manualTimeField%dEnd%sSec',fi,currFieldDirection)
                %if(fi>1)  
                %    disp('')
                %end
                if(currFieldAvgDur==0)
                    currFieldAvgDur=NaN;
                end
            else
                currFieldAvgDur=avgFieldDurations(currFieldNum);
            end
            
            localNormTimeInFieldPerTime(ti)=localTimeInFieldPerTime(ti)/currFieldAvgDur;    
        end
        
        
        directionSpecificStats.(currFieldDirection).expectedWidthPerField=expectedWidthPerField;
        directionSpecificStats.(currFieldDirection).minFieldDuration=minFieldDurations;
        directionSpecificStats.(currFieldDirection).isGoodLapForEachField=isGoodLapForEachField;
        directionSpecificStats.(currFieldDirection).avgFieldDurations=avgFieldDurations;
        directionSpecificStats.(currFieldDirection).avgFieldWidths=avgFieldWidths;
        directionSpecificStats.(currFieldDirection).fieldDurationsPerLap=fieldDurationsPerLap;
        directionSpecificStats.(currFieldDirection).fieldWidthsPerLap=fieldWidthsPerLap;
        
        directionSpecificStats.(currFieldDirection).fieldEntryPosPerTime=fieldEntryPosPerTime;
        directionSpecificStats.(currFieldDirection).fieldExitPosPerTime=fieldExitPosPerTime;
        directionSpecificStats.(currFieldDirection).fieldEntryTimePerTime=fieldEntryTimePerTime;
        directionSpecificStats.(currFieldDirection).fieldExitTimePerTime=fieldExitTimePerTime;
        directionSpecificStats.(currFieldDirection).fieldNumPerTime=fieldNumPerTime;
        directionSpecificStats.(currFieldDirection).lapNumberPerTime=lapNumberPerTime;
        
        if(strcmp(currFieldDirection,'leftward'))
            localSpikePhasePerTime(posInfo.speedAlongTrackPerTimeMSec>0)=NaN;
        else
            localSpikePhasePerTime(posInfo.speedAlongTrackPerTimeMSec<0)=NaN;
        end
        
        directionSpecificStats.(currFieldDirection).localSpikePhasePerTime=localSpikePhasePerTime;
        directionSpecificStats.(currFieldDirection).localTimeInFieldPerTime=localTimeInFieldPerTime;
        directionSpecificStats.(currFieldDirection).localDistInFieldPerTime=localDistInFieldPerTime;
        directionSpecificStats.(currFieldDirection).localNormTimeInFieldPerTime=localNormTimeInFieldPerTime;
        directionSpecificStats.(currFieldDirection).localNormDistInFieldPerTime=localNormDistInFieldPerTime;

        if(getLocalGaussSpikeRate)
            directionSpecificStats.(currFieldDirection).localGaussSpikeRatePerTime=localGaussSpikeRatePerTime(1:decFactor:end);

            directionSpecificStats.(currFieldDirection).gaussTimeAxis=gaussTimeAxis(1:decFactor:end);
       
            if((length(gaussTimeAxis)*2)==length(localGaussSpikeRatePerTime))
                localGaussSpikeRatePerTime=localGaussSpikeRatePerTime(1:2:end);
            end
      
        
            if(di==2)
                leftLocalGaussSpikeRatePerTime=localGaussSpikeRatePerTime;
                leftLocalGaussSpikeRatePerPosTime=interp1(gaussTimeAxis,leftLocalGaussSpikeRatePerTime,positionTimeAxis);
            else
                rightLocalGaussSpikeRatePerTime=localGaussSpikeRatePerTime;
                rightLocalGaussSpikeRatePerPosTime=interp1(gaussTimeAxis,rightLocalGaussSpikeRatePerTime,positionTimeAxis);
            end
        end
        %directionSpecificStats.(currFieldDirection).localGaussSpikeRatePerPosTime=localGaussSpikeRatePerPosTime;
       
        
        %{
        if(isCircMaze)
            circumference=max(currData.positionPerTimeStep);
            figure; 
            plot(positionTimeAxis,mod(currData.directionSpecificStats.(currFieldDirection).fieldEntryPosPerTime+currData.medianPos,circumference),'b','LineWidth',5)
            hold on
            plot(positionTimeAxis,mod(currData.directionSpecificStats.(currFieldDirection).fieldExitPosPerTime+currData.medianPos,circumference),'r','LineWidth',5)
            
            hold on; plot(positionTimeAxis, currData.positionPerTimeStep)
            hold on; plot(positionTimeAxis,localDistInFieldPerTime)
            hold on; plot(positionTimeAxis,currData.directionSpecificStats.(currFieldDirection).fieldNumPerTime,'k','LineWidth',5)
             hold on; plot(positionTimeAxis,fieldEntryPosPerTime,'m','LineWidth',5)
            disp('')
        end
        %}
        
        %{
        figure
        plotQuickViewSpacetime
     
        maxFig
        setFigFontTo(24)
        saveas(gcf,fullfile(spaceTimePlotSaveDir,sprintf('worldLineVsSpikePhase_%s_%s.tif',fileBaseName, currFieldDirection)))
        %pause(1)
           disp('')
        %}
           
        %{
        figure
        scatter(localDistInFieldPerTime,localTimeInFieldPerTime,20,localSpikePhasePerTime,'filled')
         %scatter(localNormDistInFieldPerTime,localNormTimeInFieldPerTime,20,localSpikePhasePerTime,'filled')
        hold on
        colormap(jet)
        colorbar
        xlabel('Distance in field (m)')
        ylabel('Time in field (sec)')
        title(removeUnderscores(getFileNameFromPath(currDataFilePath)))
        
        %avgFieldWidths=data.directionSpecificStats.(currFieldDirection).avgFieldWidths;
        %avgFieldDurations=data.directionSpecificStats.(currFieldDirection).avgFieldDurations;
        
        for i=1:length(avgFieldWidths)
            plot([abs(avgFieldWidths(i)) abs(avgFieldWidths(i))],ylim,'k--','LineWidth',3)
            plot(xlim,[avgFieldDurations(i) avgFieldDurations(i)],'k--','LineWidth',3)
        end
        
        %xlim([0 1])
        %ylim([0 1.5])
        %if(~isnan(avgFieldDurations(i)))
         %ylim([0 avgFieldDurations(i)*1.2])
        %end
        daspect([1 1 1])
        drawnow
      %}
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %plot details of field, space, time detections
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %{
        figure; plot(currData.positionTimeAxis, currData.directionSpecificStats.(currFieldDirection).localDistInFieldPerTime)
    hold on; plot(currData.positionTimeAxis,currData.directionSpecificStats.(currFieldDirection).localTimeInFieldPerTime)
    hold on;plot(currData.positionTimeAxis,fieldNumPerTime,'m-','LineWidth',2)
     plot(currData.positionTimeAxis,currData.positionPerTimeStep,'r-','LineWidth',2)
     %yyaxis right;plot(gaussTimeAxis,localGaussSpikeRatePerTime,'k-','LineWidth',2)

 fieldStartBoundDirStr=sprintf('manualFieldStartsM%s',currFieldDirection);
 fieldEndBoundDirStr=sprintf('manualFieldEndsM%s',currFieldDirection);
 firstFieldStart=currData.(fieldStartBoundDirStr)(1);
 firstFieldEnd=currData.(fieldEndBoundDirStr)(1);
plot(xlim,[firstFieldStart firstFieldStart],'g--','LineWidth',3)
plot(xlim,[firstFieldEnd firstFieldEnd],'r--','LineWidth',3)
 
plotOriginalSpikeTimes(currData)
plotLeftAndRightSpikeTimes(currData)

%try
if(length(currData.gaussTimeAxis)==length(currData.leftLocalGaussSpikeRatePerTime))
     plot(currData.gaussTimeAxis,currData.leftLocalGaussSpikeRatePerTime(1:1:end),'b-','LineWidth',2)
else
    plot(currData.gaussTimeAxis,currData.leftLocalGaussSpikeRatePerTime(1:2:end),'b-','LineWidth',2)
end
%end
%try
if(length(currData.gaussTimeAxis)==length(currData.rightLocalGaussSpikeRatePerTime))
     plot(currData.gaussTimeAxis,currData.rightLocalGaussSpikeRatePerTime(1:1:end),'m-','LineWidth',2)
else
    plot(currData.gaussTimeAxis,currData.rightLocalGaussSpikeRatePerTime(1:2:end),'m-','LineWidth',2)
end
 
%end
%plot(gaussTimeAxis,localGaussSpikeRatePerTime,'k-','LineWidth',2)

plot(currData.unitInfo.positionPerTime.TimeStamps,currData.unitInfo.positionPerTime.OneDLocation,'k--')

figure;
testPhase=currData.directionSpecificStats.(currFieldDirection).localSpikePhasePerTime;
allPhasesPerPosH=plot(currData.positionPerTimeStep,testPhase,'y.')
hold on
testPhase(isnan(fieldNumPerTime))=NaN;
allPhasesPerPosH.Color(4)=0.1;

testPos=currData.positionPerTimeStep;


%testPos(isnan(fieldNumPerTime))=NaN;

%plot(mod(currData.positionPerTimeStep+currData.medianPos,max(currData.positionPerTimeStep)),testPhase,'k.')
plot(currData.positionPerTimeStep,testPhase,'k.')
plot(currData.positionPerTimeStep,fieldNumPerTime,'b.')

figure
plot(currData.leftSpikePositions,currData.leftSpikePhases,'b.')
hold on
plot(currData.rightSpikePositions,currData.rightSpikePhases,'m.')
%figure
%plot(currData.positionTimeAxis,testPos,'k.')


%for s=1:length(currData.leftSpikeTimes(:))
%    plot([currData.leftSpikeTimes(s) currData.leftSpikeTimes(s)],ylim,'k-','LineWidth',1)
%end
        %close all
    
      autoArrangeFigures(1,2,1)
      
    disp(i)
    disp('')
   close all
   
      %}
      close all
      end %direction loop
      
    %end%try
 
    if(getLocalGaussSpikeRate)
        if(~exist('leftLocalGaussSpikeRatePerPosTime','var'))
            leftLocalGaussSpikeRatePerTime=NaN;
            leftLocalGaussSpikeRatePerPosTime=NaN;
        end
        if(~exist('rightLocalGaussSpikeRatePerPosTime','var'))
            rightLocalGaussSpikeRatePerPosTime=NaN;
            rightLocalGaussSpikeRatePerTime=NaN;
        end
        %save(currDataFilePath,'leftLocalGaussSpikeRatePerPosTime','rightLocalGaussSpikeRatePerPosTime','directionSpecificStats','rightLocalGaussSpikeRatePerTime','leftLocalGaussSpikeRatePerTime','gaussTimeAxis','-append')
        save(currDataFilePath,'directionSpecificStats','-append')

    else
         save(currDataFilePath,'directionSpecificStats','-append')
    end
    
    directionSpecificStats=[];
end %file loop

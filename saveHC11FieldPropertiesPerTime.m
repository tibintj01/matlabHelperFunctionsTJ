close all; clear all; clc
dataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles';
showPlots=0;
debugTimeFilePath='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles/Achilles_11012013_6-12Theta_Unit1024Info.mat';
debugTimeFilePath='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles/Achilles_11012013_6-12Theta_Unit1006Info.mat';
debugTimeFilePath='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles/Achilles_11012013_6-12Theta_Unit1024Info.mat';
debugTimeFilePath='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles/Achilles_11012013_6-12Theta_Unit1038Info.mat';
debugTimeFilePath='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles/Achilles_11012013_6-12Theta_Unit1114Info.mat';
debugTimeFilePath='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles/Cicero_09102014_6-12Theta_Unit929Info.mat';
debugTimeFilePath='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles/Achilles_11012013_6-12Theta_Unit1006Info.mat';
debugTimeFilePath='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles/Achilles_11012013_6-12Theta_Unit102Info.mat';
debugTimeFilePath='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles/Achilles_10252013_6-12Theta_Unit1115Info.mat';

filePaths=getFilePathsRegex(dataDir,'*mat');
totalNumFields=0;

resetGaussSpikeWidth=1;
  %kernelWidthSec=0.1; %sec
    kernelWidthSec=0.01; %sec
useManualTimeFieldBoundForNorm=1;
for i=1:length(filePaths)
%for i=1:1
    
    %try
        
    i/length(filePaths)
    currDataFilePath=filePaths{i};
    
    currFileName=getFileNameFromPath(currDataFilePath);
    fileBaseName=currFileName(1:(end-4));
    if(~strcmp(debugTimeFilePath,currDataFilePath))
       %continue
    end
    
    %skip badly detected field
     if(contains(currFileName,'Gatsby_08282013_6-12Theta_Unit210Info') && fi==3)
                    avgFieldDurations(fi)=NaN;
                end
    currData=load(currDataFilePath);
    
   
    %{
    if(max(currData.numFieldsPerDir)<2)
        continue
    end
    if(isfield(currData,'directionSpecificStats'))
        %continue
    end
    %}
    
    entryTimePerDirPerFieldPerLap=currData.entryTimePerDirPerFieldPerLap;
    exitTimePerDirPerFieldPerLap=currData.exitTimePerDirPerFieldPerLap;
    
    entryPosPerDirPerFieldPerLap=currData.entryPosPerDirPerFieldPerLap;
    exitPosPerDirPerFieldPerLap=currData.exitPosPerDirPerFieldPerLap;
    
    entryTimeIdxesPerDirPerFieldPerLap=currData.entryTimeIdxesPerDirPerFieldPerLap;
    exitTimeIdxesPerDirPerFieldPerLap=currData.exitTimeIdxesPerDirPerFieldPerLap;
    
    
  
      directionSpecificStats=[];
    for di=1:2
        if(showPlots)
            close all
              figure
          
        end
        currFieldDirection='';
        if(di==2)
            currFieldDirection='leftward';  
            if(showPlots)
               if(iscell(currData.imagePaths))
                   imshow(currData.imagePaths{1})
               else
                   imshow(currData.imagePaths)
               end
            end
        elseif(di==1)
            currFieldDirection='rightward';  
            if(showPlots)
                if(iscell(currData.imagePaths))
                   imshow(currData.imagePaths{2})
               else
                   imshow(currData.imagePaths)
                end
            end
        end
        if(showPlots)
         autoArrangeFigures(1,2,1);
        end
        if(~(isfield(currData,sprintf('manualFieldStartsM%s',currFieldDirection)) && isfield(currData,sprintf('manualFieldEndsM%s',currFieldDirection))))
            continue
        end
        
        fieldStartPosFieldName=sprintf('manualFieldStartsM%s',currFieldDirection);
        fieldEndPosFieldName=sprintf('manualFieldEndsM%s',currFieldDirection);

        startPositions=currData.(fieldStartPosFieldName);
       endPositions=currData.(fieldEndPosFieldName);
   
       cycleData=load(currData.unitInfo.cycleInfoFileName);
       cycleStartTimes=cycleData.cycleMaxTimes;
       positionTimeAxis=currData.positionTimeAxis;
       
       cycleStartTimes=cycleStartTimes(cycleStartTimes>=min(positionTimeAxis) & cycleStartTimes<=max(positionTimeAxis));
       
       if(di==1)
           allSpikeTimesAndPhases=[currData.rightSpikeTimes(:) currData.rightSpikePhases(:)];
       elseif(di==2)
            allSpikeTimesAndPhases=[currData.leftSpikeTimes(:) currData.leftSpikePhases(:)];
       end
       
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %local spike phase based on gaussian kernel peak per theta cycle
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %cycleStartTimes=
     
      
       %[localCycleMeanSpikePhasePerTime, localCycleMeanSpikeTimePerTime]=getLocalCycleMeanPhasePerTime(allSpikeTimesAndPhases,positionTimeAxis,cycleStartTimes);
     
       if(~resetGaussSpikeWidth && isfield(currData,directionSpecificStats) && isfield(currData.directionSpecificStats,currFieldDirection) && isfield(currData.directionSpecificStats.(currFieldDirection),'localSpikePhasePerTime'))
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
          [gaussTimeAxis,rightLocalGaussSpikeRatePerTime]=localGaussRateSpikeTrain(currData.rightSpikeTimes(:)-min(positionTimeAxis), kernelWidthSec,0,max(positionTimeAxis)-min(positionTimeAxis));
             gaussTimeAxis=gaussTimeAxis+min(positionTimeAxis);

             [gaussTimeAxis,leftLocalGaussSpikeRatePerTime]=localGaussRateSpikeTrain(currData.leftSpikeTimes(:)-min(positionTimeAxis), kernelWidthSec,0,max(positionTimeAxis)-min(positionTimeAxis));
             gaussTimeAxis=gaussTimeAxis+min(positionTimeAxis);

             if(di==1)
                 localGaussSpikeRatePerTime=rightLocalGaussSpikeRatePerTime;
             else
                 localGaussSpikeRatePerTime=leftLocalGaussSpikeRatePerTime;
             end

           localSpikePhasePerTime=getPeakPhasePerThetaCyclePerTime(gaussTimeAxis,localGaussSpikeRatePerTime,cycleStartTimes,positionTimeAxis);
             decFactor=2;
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
           fieldWidthsPerLap(fi,:)=currFieldWidths;
           
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
       if(max(abs(fieldWidthsPerLap(fi,:)))==0)
           continue
       end
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %create flag array indicating laps where field length is within 15 percent of expected and
       %field duration is greater than expected length times max speed 
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       if(di==1)
           expectedWidthPerField=currData.(fieldEndPosFieldName)-currData.(fieldStartPosFieldName);
       elseif(di==2)
           expectedWidthPerField=-(currData.(fieldStartPosFieldName)-currData.(fieldEndPosFieldName));
       end
       
       %{
       if(expectedWidthPerField<0)
           disp('') 
     end
       %}
       
       
       expectedWidthPerField=expectedWidthPerField(:);
       expectedWidthPerField=repmat(expectedWidthPerField,1,size(fieldWidthsPerLap,2));
       
       removeWidthRows=[];
       for wr=1:size(expectedWidthPerField,1)
           if(sum(isnan(expectedWidthPerField(wr,:)))==size(expectedWidthPerField,2))
               removeWidthRows=[removeWidthRows wr];
           end
       end
       expectedWidthPerField(removeWidthRows,:)=[];
       
       
       %thresholdWidthDeviations=expectedWidthPerField*0.15;
        thresholdWidthDeviations=0.05; %within 5 cm of expected width
        thresholdWidthDeviations=0.1; %within 10 cm of expected width
       
         %maxSpeed=2;
       maxSpeed=1.3;
       minFieldDurations=abs(expectedWidthPerField)/maxSpeed;
       
       goodFieldWidthLaps=(abs(fieldWidthsPerLap-expectedWidthPerField)<=thresholdWidthDeviations);
       goodFieldDurationLaps=fieldDurationsPerLap>minFieldDurations;
       
        isGoodLapForEachField=goodFieldWidthLaps & goodFieldDurationLaps;
        
        if(max(isGoodLapForEachField)==0)
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
                
                if(contains(currFileName,'Achilles_11012013_6-12Theta_Unit102Info') && fi==1)
                    avgFieldDurations(fi)=typicalDur*0.65;
                end
                
                 if(contains(currFileName,'Achilles_11012013_6-12Theta_Unit406Info') && fi==1)
                    avgFieldDurations(fi)=typicalDur*0.8;
                 end
                 
                 if(contains(currFileName,'Achilles_11012013_6-12Theta_Unit407Info') && fi==1)
                    avgFieldDurations(fi)=typicalDur*0.8;
                 end
                
                 if(contains(currFileName,'Achilles_11012013_6-12Theta_Unit409Info') && fi==1)
                    avgFieldDurations(fi)=typicalDur*0.65;
                 end
                if(contains(currFileName,'Achilles_11012013_6-12Theta_Unit415Info') && fi==1)
                    avgFieldDurations(fi)=typicalDur*1.4;
                end
                
                if(contains(currFileName,'Achilles_11012013_6-12Theta_Unit812Info') && fi==1)
                    avgFieldDurations(fi)=typicalDur*0.23;
                end
                if(contains(currFileName,'Achilles_11012013_6-12Theta_Unit914Info') && fi==1)
                    avgFieldDurations(fi)=typicalDur*2;
                end
                 if(contains(currFileName,'Buddy_06272013_6-12Theta_Unit215Info') && di==1)
                    avgFieldDurations(fi)=typicalDur*0.9;
                 end
                if(contains(currFileName,'Buddy_06272013_6-12Theta_Unit219Info') && di==2)
                    avgFieldDurations(fi)=typicalDur*0.9;
                end
                 if(contains(currFileName,'Buddy_06272013_6-12Theta_Unit605Info') && di==2)
                    avgFieldDurations(fi)=typicalDur*0.85;
                 end
                if(contains(currFileName,'Cicero_09102014_6-12Theta_Unit615Info') && di==2)
                    avgFieldDurations(fi)=typicalDur*0.85;
                end
                if(contains(currFileName,'Cicero_09102014_6-12Theta_Unit817Info') && di==2)
                    avgFieldDurations(fi)=typicalDur*0.85;
                end
                if(contains(currFileName,'Cicero_09102014_6-12Theta_Unit105Info') && fi==2)
                    avgFieldDurations(fi)=NaN;
                end
                 if(contains(currFileName,'Cicero_09102014_6-12Theta_Unit607Info') && fi==1)
                    avgFieldDurations(fi)=NaN;
                 end
                if(contains(currFileName,'Gatsby_08282013_6-12Theta_Unit206Info') && fi==1)
                    avgFieldDurations(fi)=NaN;
                end
                if(contains(currFileName,'Gatsby_08282013_6-12Theta_Unit719Info') && di==2)
                    avgFieldDurations(fi)=typicalDur*0.85;
                end
                if(contains(currFileName,'Gatsby_08282013_6-12Theta_Unit721Info') && di==2)
                    avgFieldDurations(fi)=typicalDur*0.85;
                end
                if(contains(currFileName,'Gatsby_08282013_6-12Theta_Unit210Info') && fi==3)
                    avgFieldDurations(fi)=NaN;
                end
                 if(contains(currFileName,'Gatsby_08282013_6-12Theta_Unit608Info') && fi==3)
                    avgFieldDurations(fi)=NaN;
                end
                %avgFieldDurations(fi)=nanmean((goodFieldDurs));
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
        localDistInFieldPerTime=NaN(size(currData.positionPerTimeStep));
        localTimeInFieldPerTime=NaN(size(currData.positionTimeAxis));
        
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
            
            if(isCircMaze)
                circumference=max(currData.positionPerTimeStep);
                %adjustedEntryPos=mod(fieldEntryPosPerTime(ti)-currData.medianPos,circumference);
                adjustedEntryPos=fieldEntryPosPerTime(ti);
                %shiftPos=mod(currData.positionPerTimeStep(ti)+currData.medianPos,max(currData.positionPerTimeStep));
                localDistInFieldPerTime(ti)=abs(circDiffGivenCircumference(currData.positionPerTimeStep(ti),adjustedEntryPos,max(currData.positionPerTimeStep)));
                %localDistInFieldPerTime(ti)=abs(circDiffGivenCircumference(shiftPos,adjustedEntryPos,max(currData.positionPerTimeStep)));

            else
                localDistInFieldPerTime(ti)=abs(currData.positionPerTimeStep(ti)-fieldEntryPosPerTime(ti));
            end
            
             localTimeInFieldPerTime(ti)=abs(currData.positionTimeAxis(ti)-fieldEntryTimePerTime(ti));
            
            currFieldAvgWidth=avgFieldWidths(currFieldNum);
            localNormDistInFieldPerTime(ti)=localDistInFieldPerTime(ti)/abs(currFieldAvgWidth);
            
            if(useManualTimeFieldBoundForNorm)
                %currFieldNum reset for each time point
                currFieldManualDuration=currData.(sprintf('manualTimeField%dEnd%sSec',currFieldNum,currFieldDirection));
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
            localSpikePhasePerTime(currData.speedPerTimeStep>0)=NaN;
        else
            localSpikePhasePerTime(currData.speedPerTimeStep<0)=NaN;
        end
        
        directionSpecificStats.(currFieldDirection).localSpikePhasePerTime=localSpikePhasePerTime;
        
        directionSpecificStats.(currFieldDirection).localTimeInFieldPerTime=localTimeInFieldPerTime;
        directionSpecificStats.(currFieldDirection).localDistInFieldPerTime=localDistInFieldPerTime;
        
        directionSpecificStats.(currFieldDirection).localNormTimeInFieldPerTime=localNormTimeInFieldPerTime;
        directionSpecificStats.(currFieldDirection).localNormDistInFieldPerTime=localNormDistInFieldPerTime;

        directionSpecificStats.(currFieldDirection).localGaussSpikeRatePerTime=localGaussSpikeRatePerTime(1:decFactor:end);
        directionSpecificStats.(currFieldDirection).gaussTimeAxis=gaussTimeAxis(1:decFactor:end);

        if((length(gaussTimeAxis)*2)==length(leftLocalGaussSpikeRatePerTime))
            leftLocalGaussSpikeRatePerTime=leftLocalGaussSpikeRatePerTime(1:2:end);
        end
        if((length(gaussTimeAxis)*2)==length(rightLocalGaussSpikeRatePerTime))
            rightLocalGaussSpikeRatePerTime=rightLocalGaussSpikeRatePerTime(1:2:end);
        end
        leftLocalGaussSpikeRatePerPosTime=interp1(gaussTimeAxis,leftLocalGaussSpikeRatePerTime,currData.positionTimeAxis);
        rightLocalGaussSpikeRatePerPosTime=interp1(gaussTimeAxis,rightLocalGaussSpikeRatePerTime,currData.positionTimeAxis);
        %directionSpecificStats.(currFieldDirection).localGaussSpikeRatePerPosTime=localGaussSpikeRatePerPosTime;
       
        save(currDataFilePath,'leftLocalGaussSpikeRatePerPosTime','rightLocalGaussSpikeRatePerPosTime','directionSpecificStats','rightLocalGaussSpikeRatePerTime','leftLocalGaussSpikeRatePerTime','gaussTimeAxis','-append')
        
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
        
        figure
        plotQuickViewSpacetime
     
        
        %pause(1)
           disp('')
           
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
  
end %file loop

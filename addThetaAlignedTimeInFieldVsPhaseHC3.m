close all; clear all; clc
processedDataDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData';
unitDataDir=fullfile(processedDataDir,'unitSpikeInfoStructs');

savePhaseVsThetaCyclesElapsedDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3Images/phaseVsTimeInFieldInCycles';

touchDir(savePhaseVsThetaCyclesElapsedDir)

filePaths=getFilePathsRegex(unitDataDir,'*mat');

firstSpikeAsStartOfTimeField=1;
firstSpikeAsStartOfTimeField=0;

colorCodeForSpeed=0;
colorCodeForPhase=1;
showPlots=0;
%showPlots=1;

showRaster=0;

minNumSpikesPerTraversal=3;
justFirstSpikes=1;
justFirstSpikes=0;

useRefThetaCh=0;
useRefThetaCh=1;

maxTimeDisp=1.5;
maxTimeDisp=3;
%autoShiftOffset=1;

cmapName=magma(100);
%cmapName=plasma(100);
speedColors=cmapName;

%minAbsSpeedInLap=0.05;
minAbsSpeedInLap=0.01;
minAbsSpeedInLap=0;
maxFirstCycleInField=5;
maxFirstCycleInField=3;
%minLapNumSpikes=3;
%minLapNumSpikes=5;
minLapNumSpikes=3;
eliminateBackwards=1;


phaseColors=hsv(360);
%phaseColors=phasemap(360);

startFile=1;
%{
startFile=48;
startFile=84;
%}
%startFile=75;
for fi=startFile:length(filePaths)
    %for fi=76:length(filePaths)
       fi
    currFilePath=filePaths{fi};
    currFileName=getFileNameFromPath(currFilePath);
    fileBaseName=currFileName(1:(end-4));
    thetaBasedTimeVsPhaseInfo=[];

    %{
    if(~strcmp(currFilePath,'/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData/unitSpikeInfoStructs/ec014.468_Unit102Info.mat'))
        continue
    end
   %}
    
    %{
    [currSessionNum,currSessName]=getSessionNum(fileBaseName);
    
    if(isnan(currSessName))
        continue
    end
    
    
    if(~contains(currSessName,'atsby'))
        %continue
    end
    
     [leftwardRefUnitData,rightwardRefUnitData]=getRefUnitData(currSessName);
    %}
    
    dataStruct=load(currFilePath);
    data=dataStruct.unitInfo;
    
    allLapInFieldSpikeTimes=[];
    allLapInFieldSpikeTimesInExp=[];
    allLapInFieldSpikeTraversalAvgSpeeds=[];
     allLapInFieldSpikeTraversalSpeedRanges=[];
    
    
               allLapInFieldSpikePhases=[];
               allLapSpikeCycleIDsFromFieldEntry=[];
               allLapInFieldCycleMinTimes=[];
               allLapInFieldSpikeSpeeds=[];
               allLapInFieldSpikeSpeedsSigned=[];
               allLapInFieldMeanCycleDurs=[];
               allLapInFieldSpikesInLapWithLateStart=[];
               
               allLapInFieldMeanPhasePerTraversal=[];
               allLapInFieldMeanSpeedPerTraversal=[];
               allLapInFieldSpeedRangePerTraversal=[];
               

   if(~isfield(dataStruct,'entryTimePerDirPerFieldPerLap'))
       continue
   end
    entryTimePerDirPerFieldPerLap=dataStruct.entryTimePerDirPerFieldPerLap;
    exitTimePerDirPerFieldPerLap=dataStruct.exitTimePerDirPerFieldPerLap;
    
    for di=1:2
            %figure(1)
            
            allLapInFieldSpikeTimes=[];
            allLapInFieldSpikeTimesInExp=[];
            allLapInFieldSpikeTraversalAvgSpeeds=[];
            allLapInFieldSpikeTraversalSpeedRanges=[];
            allLapInFieldSpikesInLapWithLateStart=[];
               allLapInFieldSpikePhases=[];
               allLapSpikeCycleIDsFromFieldEntry=[];
               allLapInFieldCycleMinTimes=[];
               allLapInFieldSpikeSpeeds=[];
               
               allLapInFieldMeanPhasePerTraversal=[];
               allLapInFieldMeanSpeedPerTraversal=[];
               allLapInFieldSpeedRangePerTraversal=[];
               
               
        if(di==2)
            currFieldDirection='leftward'; 
            goingBackwardIdxes=data.speedPerTimeStep>0;
            %gaussSpikeRatePerTime=data.leftLocalGaussSpikeRatePerPosTime;
            %refUnitData=leftwardRefUnitData;
            

             spikePhases=data.leftSpikePhases;
            %spikeTimes=data.leftSpikeTimes;
            spikeTimes=data.spikeTimes(data.spikeTimeDirectionAssignment==di);
            spikePositions=data.leftSpikePositions;
            
        elseif(di==1)
            currFieldDirection='rightward'; 
            goingBackwardIdxes=data.speedPerTimeStep<0;
            %refUnitData=rightwardRefUnitData;
                      
            %gaussSpikeRatePerTime=data.rightLocalGaussSpikeRatePerPosTime;
            
            spikePhases=data.rightSpikePhases;
            %spikeTimes=data.rightSpikeTimes;
            spikeTimes=data.spikeTimes(data.spikeTimeDirectionAssignment==di);
            spikePositions=data.rightSpikePositions;
        end
       
    spikePhasesPreShift=spikePhases;
    %{
     if(~useRefThetaCh || isnan(data.refThetaCh))
        spikePhases=manualApproxShiftDeg(spikePhasesPreShift,currSessName);
     else 
        spikePhases=manualApproxShiftRefDeg(spikePhasesPreShift,currSessName);
     end 
     %}
        posInfo=load(data.positionInfoForSessionSavePath);
        positionTimeAxis=posInfo.positionTimeAxisSec;
        
        spikeSpeeds=abs(interp1(positionTimeAxis,dataStruct.filledSignedSpeedPerTime,spikeTimes));
        signedSpikeSpeeds=interp1(positionTimeAxis,dataStruct.filledSignedSpeedPerTime,spikeTimes);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %load space time coordinates in field per time
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %{
        if(~isfield(data.directionSpecificStats,currFieldDirection))
            close all
            continue
        end

        numFields=size(data.directionSpecificStats.(currFieldDirection).expectedWidthPerField,1);
        %}
        
        numFields=length(dataStruct.(sprintf('%sFieldStartEndM',currFieldDirection)))/2;
        
        if(~isfield(dataStruct.lapStartTimesPerDir,currFieldDirection))
            continue
        end
        lapStartTimes=dataStruct.lapStartTimesPerDir.(currFieldDirection);
        lapEndTimes=dataStruct.lapStopTimesPerDir.(currFieldDirection);
        
        refCycleInfo=load(dataStruct.refChCycleInfoPath);
        refCycleMaxTimes=refCycleInfo.cycleMaxTimes;
        
        numLaps=length(lapStartTimes);
        
        %positionTimeAxis=data.positionTimeAxis;
        
        refCycleMaxTimes=refCycleMaxTimes(refCycleMaxTimes>=positionTimeAxis(1) & refCycleMaxTimes<=positionTimeAxis(end));
        
        approxTimeStep=median(diff(positionTimeAxis));
        maxTimeSince=5;%seconds
        maxTimeSince=1.5;%seconds
         maxTimeSince=0.75;%seconds
         maxDistSince=1; %cm
         
        timeSinceAxis=(approxTimeStep/2):approxTimeStep:maxTimeSince;
        numTimePts=length(timeSinceAxis);
        
        %gaussRateProfile=NaN(numLaps,numTimePts);
        
        currDirBounds=dataStruct.(sprintf('%sFieldStartEndM',currFieldDirection));
        
        currFieldStartPosPerField=NaN(numFields,1);
        currFieldEndPosPerField=NaN(numFields,1);
        
        fieldCount=0;
        for bi=1:2:length(currDirBounds)
            fieldCount=fieldCount+1;
            currFieldStartPosPerField(fieldCount)=currDirBounds(bi);
            currFieldEndPosPerField(fieldCount)=currDirBounds(bi+1);
        end
          
        %currFieldStartPosPerField=data.(sprintf('%sFieldStartEndM',currFieldDirection));
        %currFieldEndPosPerField=data.(sprintf('manualFieldEndsM%s',currFieldDirection));

        currDirectionTimeBounds=dataStruct.cycleTimeFieldBoundPerDirPerField.(currFieldDirection)-1; %-1 becaues manually entered before changing first theta cycle from 1 to 0
        for fii=1:numFields
            
           saveImageFilePath=fullfile(savePhaseVsThetaCyclesElapsedDir,sprintf('sinceFieldEntrySpikePhaseVsTheta_%s_%s_Field%d.tif',fileBaseName,currFieldDirection,fii));

           
           
            
           if(exist(saveImageFilePath,'file'))
               % continue
            end
            
            allLapInFieldSpikeTimes=[];
            allLapInFieldSpikeTimesInExp=[];
            allLapInFieldSpikeTraversalAvgSpeeds=[];
            allLapInFieldSpikeTraversalSpeedRanges=[];
            allLapInFieldSpikesInLapWithLateStart=[];
            allLapInFieldSpikeLapNums=[];
               allLapInFieldSpikePhases=[];
               allLapSpikeCycleIDsFromFieldEntry=[];
               allLapInFieldCycleMinTimes=[];
               allLapInFieldSpikeSpeeds=[];
               allLapInFieldSpikeSpeedsSigned=[];
               allLapInFieldSpikesInLapWithStop=[];
               allLapInFieldSpikeDists=[];
               allLapInFieldSpikeDistsFieldFrac=[];
               allLapInFieldNumCycles=[];
               
               currFieldStr=sprintf('field%d',fii);
               
               if(showPlots)
                   if(showRaster)
                    fRaster=figure(20);
                   end
                    fPhaseVsCycleNum=figure(21);
                    %autoArrangeFigures()
               end
               %{
            if(fii>length(data.directionSpecificStats.(currFieldDirection).avgFieldDurations))
                continue
            end
            %}
               
               %{
            try
                fieldDurLim=1.5*data.directionSpecificStats.(currFieldDirection).avgFieldDurations(fii);
                fieldEntryTimePerTime=data.directionSpecificStats.(currFieldDirection).fieldEntryTimePerTime(fii,:);
                fieldExitTimePerTime=data.directionSpecificStats.(currFieldDirection).fieldExitTimePerTime(fii,:);
                fieldWidth=abs(data.directionSpecificStats.(currFieldDirection).avgFieldWidths(fii));
                fieldDur=abs(data.directionSpecificStats.(currFieldDirection).avgFieldDurations(fii));
            catch
                continue
            end
               %}
            
            %{
            if(autoShiftOffset)
                [bestShiftDeg]=getBestBySpaceShiftDeg(spikePhases,(spikePositions-min(spikePositions))/fieldWidth);
                spikePhases=mod(spikePhases-bestShiftDeg,360);
            else
            %}
                
            %end
        
            currFieldStartPos=currFieldStartPosPerField(fii);
            currFieldEndPos=currFieldEndPosPerField(fii);

	    currFieldWidthM=abs(currFieldStartPos-currFieldEndPos);
 
            numThetaCyclesInFieldPerLap=NaN(numLaps,1);
            
            allLapInFieldAllTimeVsPhasePerTraversal=cell(1);
            traversalNum=0;
            
            
            for li=1:numLaps
                
                inLapIdxes=positionTimeAxis>=lapStartTimes(li) & positionTimeAxis < lapEndTimes(li);
                %fieldEntryTimeThisLap=nanmean(fieldEntryTimePerTime(inLapIdxes));
                if(li>size(entryTimePerDirPerFieldPerLap,3))
                    continue
                end
                fieldEntryTimeThisLap=entryTimePerDirPerFieldPerLap(di,fii,li);
                fieldExitTimeThisLap=exitTimePerDirPerFieldPerLap(di,fii,li);
                
                
                inFieldTimeIdxes=positionTimeAxis>=fieldEntryTimeThisLap & positionTimeAxis<fieldExitTimeThisLap;
                speedsInThisFieldLap=posInfo.speedAlongTrackPerTimeMSec(inFieldTimeIdxes);
                meanSpeedInThisFieldLap=nanmean(speedsInThisFieldLap);
                
                
                
                %fieldExitTimeThisLap=nanmean(fieldExitTimePerTime(inLapIdxes));
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %get closest theta trough time relative to field entry time
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [~,closestRefCycleMinIdx]=min(abs(refCycleMaxTimes-fieldEntryTimeThisLap));
                if(refCycleMaxTimes(closestRefCycleMinIdx)<fieldEntryTimeThisLap)
                    closestRefCycleMinIdx=max(1,closestRefCycleMinIdx-1);
                end
                closestThetaCycleMinTime=refCycleMaxTimes(closestRefCycleMinIdx);
                
                inLapInFieldCycleMinTimeIdxes=refCycleMaxTimes>=lapStartTimes(li) &  refCycleMaxTimes< lapEndTimes(li) & refCycleMaxTimes>=fieldEntryTimeThisLap &  refCycleMaxTimes<fieldExitTimeThisLap;
                
                inFieldCycleTimes=refCycleMaxTimes(inLapInFieldCycleMinTimeIdxes);
                inFieldCycleDurations=diff(inFieldCycleTimes);
                meanInFieldCycleDuration=nanmean(inFieldCycleDurations);
                
                
                numThetaCyclesInFieldPerLap(li)=sum(inLapInFieldCycleMinTimeIdxes);
                
                inLapSpikeTimeIdxes=spikeTimes>=lapStartTimes(li) & spikeTimes < lapEndTimes(li);
                inLapInFieldSpikeTimeIdxes= inLapSpikeTimeIdxes & spikeTimes>=fieldEntryTimeThisLap & spikeTimes < fieldExitTimeThisLap;

                if(max(inLapInFieldSpikeTimeIdxes)==0)
                    continue
                end
                inLapInFieldCycleMinTimes=refCycleMaxTimes(inLapInFieldCycleMinTimeIdxes);
                minCycleBoundTime=min(inLapInFieldCycleMinTimes);
                maxCycleBoundTime=max(inLapInFieldCycleMinTimes);
                
                inLapInFieldSpikeTimesInLap=spikeTimes(inLapInFieldSpikeTimeIdxes);
                
                   inLapInFieldSpikeTimesInExp=spikeTimes(inLapInFieldSpikeTimeIdxes);
                
                inLapInFieldSpikeTimesInLap(inLapInFieldSpikeTimesInLap<minCycleBoundTime)=NaN; %can't assign spikes prior to lowest cycle bound
                %inLapInFieldSpikeTimes(inLapInFieldSpikeTimes>maxCycleBoundTime
                
               
                inLapInFieldSpikeSpeeds=spikeSpeeds(inLapInFieldSpikeTimeIdxes);
                inLapInFieldSpikeSpeedsSigned=signedSpikeSpeeds(inLapInFieldSpikeTimeIdxes);
                inLapInFieldSpikePositions=spikePositions(inLapInFieldSpikeTimeIdxes);
                inLapInFieldSpikePhases=spikePhases(inLapInFieldSpikeTimeIdxes);
                
                inLapInFieldSpikeLapNums=repelem(li,length(inLapInFieldSpikeTimesInLap));
                
                %spikeCycleIDsFromFieldEntry=NaN(size(inLapInFieldSpikeTimes));
                
                %interpolate cycle times at spike times to get cycle IDs
                if(length(inLapInFieldCycleMinTimes)>1)
                    spikeCycleIDsFromFieldEntry=floor(interp1(inLapInFieldCycleMinTimes,(1:(length(inLapInFieldCycleMinTimes)))-1,inLapInFieldSpikeTimesInLap));
                else
                    %spikeCycleIDsFromFieldEntry=1;
                    spikeCycleIDsFromFieldEntry=0;
                end
                
                %spikeCycleIDsFromFieldEntry(isnan(spikeCycleIDsFromFieldEntry))=[]; %extra NaNs.....
               
                %inLapInFieldSpikeTimes=inLapInFieldSpikeTimes-min(inLapInFieldSpikeTimes);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %zero relative to theta cycle timing so that theta is time
                %keeper across laps
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             
                inLapInFieldSpikeTimesInLap=inLapInFieldSpikeTimesInLap-closestThetaCycleMinTime;
               %inLapInFieldSpikeTimesInLap=inLapInFieldSpikeTimes-closestThetaCycleMinTime;

                inLapInFieldCycleMinTimes=inLapInFieldCycleMinTimes-closestThetaCycleMinTime;
                
                %inLapInFieldSpikeDists=abs(inLapInFieldSpikePositions-min(inLapInFieldSpikePositions));
                if(strcmp(currFieldDirection,'rightward'))
                    outOfFieldSpikeIdxes=inLapInFieldSpikePositions<currFieldStartPos | inLapInFieldSpikePositions>currFieldEndPos;
                else
                    outOfFieldSpikeIdxes=inLapInFieldSpikePositions>currFieldStartPos | inLapInFieldSpikePositions<currFieldEndPos;
                end
                
                inLapInFieldSpikeDists=abs(inLapInFieldSpikePositions-currFieldStartPos);
                %inLapInFieldSpikeDists=inLapInFieldSpikeDists*100; %m to cm
                
                %{
                inLapInFieldSpikeDists(inLapInFieldSpikeTimes>maxTimeSince)=[];
                inLapInFieldSpikeSpeeds(inLapInFieldSpikeTimes>maxTimeSince)=[];
                inLapInFieldSpikeTimes(inLapInFieldSpikeTimes>maxTimeSince)=[];
                %}
                  
                
                %plotRasterStyle(inLapInFieldSpikeTimes,li,NaN,NaN,NaN);
                %speedIdx=min(100,ceil(nanmean(inLapInFieldSpikeSpeeds*100)));
                %tickColor=speedColors(speedIdx,:);
                if(isnan(inLapInFieldSpikeSpeeds))
                    badSpeedIdxes=1;
                else
                    badSpeedIdxes=[];
                end
                
                %badSpeedIdxes=badSpeedIdxes | inLapInFieldSpikeSpeeds<0.05;
                
                inLapInFieldSpikeSpeeds(badSpeedIdxes)=[];
                inLapInFieldSpikeSpeedsSigned(badSpeedIdxes)=[];
                inLapInFieldSpikeTimesInLap(badSpeedIdxes)=[];
                inLapInFieldSpikeDists(badSpeedIdxes)=[];
                inLapInFieldSpikePhases(badSpeedIdxes)=[];
                inLapInFieldSpikeTimesInExp(badSpeedIdxes)=[];
                
                try
                    %{
                    if(~isempty(inLapInFieldSpikeSpeeds))
                     inLapInFieldSpikeSpeeds(outOfFieldSpikeIdxes)=[];
                    end

                   if(~isempty(inLapInFieldSpikeTimes))
                    inLapInFieldSpikeTimes(outOfFieldSpikeIdxes)=[];
                   end
                   if(~isempty(inLapInFieldSpikeDists))
                    inLapInFieldSpikeDists(outOfFieldSpikeIdxes)=[];
                   end
                   if(~isempty(inLapInFieldSpikePhases))
                    inLapInFieldSpikePhases(outOfFieldSpikeIdxes)=[];
                   end
                   
                   if(~isempty(spikeCycleIDsFromFieldEntry))
                    spikeCycleIDsFromFieldEntry(outOfFieldSpikeIdxes)=[];
                   end
                   %}
                    if(~isempty(inLapInFieldSpikeSpeeds))
                     inLapInFieldSpikeSpeeds(outOfFieldSpikeIdxes)=NaN;
                    end
                    
                    if(~isempty(inLapInFieldSpikeSpeedsSigned))
                        inLapInFieldSpikeSpeedsSigned(outOfFieldSpikeIdxes)=NaN;
                    end

                   if(~isempty(inLapInFieldSpikeTimesInLap))
                    inLapInFieldSpikeTimesInLap(outOfFieldSpikeIdxes)=NaN;
                   end
                   
                    if(~isempty(inLapInFieldSpikeTimesInExp))
                   inLapInFieldSpikeTimesInExp(outOfFieldSpikeIdxes)=NaN;
                    end
                   
                   if(~isempty(inLapInFieldSpikeDists))
                    inLapInFieldSpikeDists(outOfFieldSpikeIdxes)=NaN;
                   end
                   if(~isempty(inLapInFieldSpikePhases))
                    inLapInFieldSpikePhases(outOfFieldSpikeIdxes)=NaN;
                   end
                   
                   if(~isempty(spikeCycleIDsFromFieldEntry))
                    spikeCycleIDsFromFieldEntry(outOfFieldSpikeIdxes)=NaN;
                   end
                end

                if(colorCodeForSpeed)
                    speedIdxes=ceil(inLapInFieldSpikeSpeeds*100);
                    speedIdxes(speedIdxes>100)=100;
                    speedIdxes(isnan(speedIdxes))=100;
                    tickColor=speedColors(speedIdxes,:);
                elseif(colorCodeForPhase)
                    phaseIdxes=ceil(inLapInFieldSpikePhases);
                    phaseIdxes(phaseIdxes==0)=1; %ceil for exact 0
                    
                    phaseIdxes(isnan(phaseIdxes))=1; %temp fix
                    tickColor=phaseColors(phaseIdxes,:);
                end
                
                 if(isempty(inLapInFieldSpikeTimesInLap))
                    continue
                end
                %subplot(2,1,1)
                if(showPlots)
                    figure(fPhaseVsCycleNum)
                    plot(spikeCycleIDsFromFieldEntry,inLapInFieldSpikePhases,'k.','MarkerSize',15)
                    hold on
                   xlabel('Theta cycle since field entry (num)')
                      ylabel('Theta phase (deg)')

                      xlim([0 maxTimeDisp/0.1])

                    if(showRaster)
                       figure(fRaster)

                        if(li==1 || (~exist('h1','var')) || ~isvalid(h1))
                         %h1=subplot(1,10,1:9);
                         %h1=subplot(10,1,1:9);
                          h1=subplot(10,1,1:6);
                        else
                            hold(h1,'on')
                             axes(h1)
                        end
                    end
                    
                end
                
                %plotRasterStyle(inLapInFieldSpikeTimes/fieldDur,li,NaN,NaN,tickColor);
                inLapInFieldSpikePhases=inLapInFieldSpikePhases(inLapInFieldSpikeTimesInLap<=maxTimeDisp);
                inLapInFieldSpikeSpeeds=inLapInFieldSpikeSpeeds(inLapInFieldSpikeTimesInLap<=maxTimeDisp);
                inLapInFieldSpikeSpeedsSigned=inLapInFieldSpikeSpeedsSigned(inLapInFieldSpikeTimesInLap<=maxTimeDisp);
                spikeCycleIDsFromFieldEntry=spikeCycleIDsFromFieldEntry(inLapInFieldSpikeTimesInLap<=maxTimeDisp);
                inLapInFieldSpikeDists=inLapInFieldSpikeDists(inLapInFieldSpikeTimesInLap<=maxTimeDisp);
                inLapInFieldSpikeLapNums=inLapInFieldSpikeLapNums(inLapInFieldSpikeTimesInLap<=maxTimeDisp);
                
                inLapInFieldSpikeTimesInExp=inLapInFieldSpikeTimesInExp(inLapInFieldSpikeTimesInLap<=maxTimeDisp);
               
		inLapInFieldSpikeTimesInLap=inLapInFieldSpikeTimesInLap(inLapInFieldSpikeTimesInLap<=maxTimeDisp);
		 
                inLapInFieldCycleMinTimes=inLapInFieldCycleMinTimes(inLapInFieldCycleMinTimes<=maxTimeDisp);
                numCyclesInLapInField=length(inLapInFieldCycleMinTimes);
                
                if(showRaster)
                    plotRasterStyle(inLapInFieldSpikeTimesInLap,li,NaN,NaN,tickColor);
                    colormap(gca,magma)
                    hold on
                    plotRasterStyle(inLapInFieldCycleMinTimes,li,NaN,NaN,NaN);

                    cb=colorbar('NorthOutside')
                    ylabel(cb,'running speed at spike time (cm/s)')

                    %xlim([0 maxTimeSince])
                    xlim([0 maxTimeDisp])
                    caxis([0 100])
                     if(colorCodeForPhase)
                         colormap(gca,hsv)
                         %colormap(gca,phasemap)
                         if(useRefThetaCh)
                               ylabel(cb,'Ref. channel theta phase (deg)')
                           else
                             ylabel(cb,'Local theta phase (deg)')
                           end
                        caxis([0 360])
                     end

                    end
                %drawnow
                
                 %subplot(2,1,2)
                
              
                 %subplot(2,10,11:19)
                 %{
                 if(li==1 || ~isvalid(h2))
                    h2=subplot(2,10,11:19);
                else
                    hold(h2,'on')
                     axes(h2)
                 end
                 
                 xlabel('Theta cycle since field entry (no.)')
                ylabel('Lap number')
                 %}
                
                 %{
                plotRasterStyle(inLapInFieldSpikeDists/fieldWidth,li,NaN,NaN,tickColor);
                hold on
                colormap(gca,cmapName)
                cb2=colorbar;
                ylabel(cb2,'running speed at spike time (cm/s)')
                
                xlim([0 1])
                caxis([0 100])
                if(colorCodeForPhase)
                     colormap(gca,hsv)
                     %colormap(gca,magma)
                       %colormap(gca,phasemap)
                       if(useRefThetaCh)
                           ylabel(cb2,'Ref. theta phase (deg)')
                       else
                         ylabel(cb2,'Local theta phase (deg)')
                       end
                    caxis([0 360])
                end
                
                
                %}
                  %subplot(1,10,10)
                  %subplot(10,1,10)
                  
                  if(showRaster)
                      
                      subplot(10,1,[7 8 9 10])
                    %plot(inLapInFieldSpikeTimes/fieldDur,inLapInFieldSpikePhases,'k.','MarkerSize',10)

                   plot(inLapInFieldSpikeTimesInLap,inLapInFieldSpikePhases,'k.','MarkerSize',10)


                    hold on
                    %uberTitle(removeUnderscores(currSessName))
                    %if(li==1)
                    xlim([0 maxTimeDisp])
                     ylim([0 360])

                        %xlabel('Time (frac field)')
                        xlabel('Time since field entry (sec)')
                      ylabel('Theta phase (deg)')
                        %yticklabels({})
                        %title('Time vs phase')
                    %end
                
                  end
                  
                  %if(min(spikeCycleIDsFromFieldEntry)>maxFirstCycleInField)
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  %mark this lap if it has an unusual feature like
                  %backwards movement
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  %if(length(spikeCycleIDsFromFieldEntry)<minLapNumSpikes || nanmean(spikeCycleIDsFromFieldEntry(1:minLapNumSpikes))>maxFirstCycleInField)
                  if(length(spikeCycleIDsFromFieldEntry)<minLapNumSpikes || max(spikeCycleIDsFromFieldEntry(1:minLapNumSpikes))>maxFirstCycleInField)
                      inLapInFieldSpikesInLapWithLateStart=ones(length(spikeCycleIDsFromFieldEntry),1);
                  else
                      inLapInFieldSpikesInLapWithLateStart=zeros(length(spikeCycleIDsFromFieldEntry),1);
                  end
                  
                  if(eliminateBackwards)
                      if(strcmp(currFieldDirection,'rightward'))
                            
                          if(min(inLapInFieldSpikeSpeedsSigned)<0)
                              inLapInFieldSpikesInLapWithStop=ones(length(inLapInFieldSpikeSpeeds),1);
                          else
                              inLapInFieldSpikesInLapWithStop=zeros(length(inLapInFieldSpikeSpeeds),1);
                          end
                      else
                          if(max(inLapInFieldSpikeSpeedsSigned)>0)
                              inLapInFieldSpikesInLapWithStop=ones(length(inLapInFieldSpikeSpeeds),1);
                          else
                              inLapInFieldSpikesInLapWithStop=zeros(length(inLapInFieldSpikeSpeeds),1);
                          end
                      end
                  end
                  
                  if(firstSpikeAsStartOfTimeField)
                     spikeCycleIDsFromFieldEntry=spikeCycleIDsFromFieldEntry-min(spikeCycleIDsFromFieldEntry)+1;
                  end
                      
                      %{
                  if(min(inLapInFieldSpikeSpeeds)<minAbsSpeedInLap)
                          inLapInFieldSpikesInLapWithStop=ones(length(inLapInFieldSpikeSpeeds),1);
                      else
                          inLapInFieldSpikesInLapWithStop=zeros(length(inLapInFieldSpikeSpeeds),1);
                      end
                      %}
                  %speedRangeInThisFieldLap=range(speedsInThisFieldLap);
                   speedRangeInThisFieldLap=nanstd(speedsInThisFieldLap);
                    %speedRangeInThisFieldLap=getSEMacrossRows(speedsInThisFieldLap(:));
                    
                  
                  inLapInFieldSpikeTraversalAvgSpeeds=repelem(meanSpeedInThisFieldLap,length(inLapInFieldSpikeTimesInLap));
                  inLapInFieldSpikeTraversalSpeedRanges=repelem(speedRangeInThisFieldLap,length(inLapInFieldSpikeTimesInLap));
                  
                  if(length(inLapInFieldSpikePhases(:))>=minNumSpikesPerTraversal)
                    thisLapInFieldMeanPhase=circMeanDeg(inLapInFieldSpikePhases(:));
                  else
                      thisLapInFieldMeanPhase=NaN;
                  end
                  
                  
                    
                    
                    
                    %currTimeBound=currDirectionTimeBounds(fii)*0.8*0.9; %actual bound for log 0 is 0.7 what it was
                    %currTimeBound=currDirectionTimeBounds(fii)*0.75*0.95; %actual bound for log 0 is 0.7 what it was
                    %currTimeBound=currDirectionTimeBounds(fii)*0.75; %actual bound for log 0 is 0.7 what it was
                    %currTimeBound=currDirectionTimeBounds(fii)*0.8; %actual bound for log 0 is 0.7 what it was

                  %currTimeBound=currDirectionTimeBounds(fii)*0.75*0.9; %actual bound for log 0 is 0.7 what it was
                  %currTimeBound=currDirectionTimeBounds(fii)*0.75*0.9; %actual bound for log 0 is 0.7 what it was
                  %currTimeBound=currDirectionTimeBounds(fii)*0.9; %actual bound for log 0 is 0.7 what it was
                 % currTimeBound=currDirectionTimeBounds(fii); %actual bound for log 0 is 0.7 what it was
                    currTimeBound=currDirectionTimeBounds(fii)*0.8

                %currTimeBound=currDirectionTimeBounds(fii)*0.8*0.9; %actual bound for log 0 is 0.7 what it was

                    inLapInFieldSpikePhasesTimeFracs=spikeCycleIDsFromFieldEntry/currTimeBound;
                
                     traversalNum=traversalNum+1;
                  allLapInFieldAllTimeVsPhasePerTraversal{traversalNum}=[inLapInFieldSpikePhasesTimeFracs(:) inLapInFieldSpikePhases(:)];
                  
                 
                    
                  
                  allLapInFieldMeanPhasePerTraversal=[allLapInFieldMeanPhasePerTraversal;thisLapInFieldMeanPhase];
                  
                  allLapInFieldMeanSpeedPerTraversal=[allLapInFieldMeanSpeedPerTraversal;abs(meanSpeedInThisFieldLap)];
                  allLapInFieldSpeedRangePerTraversal=[allLapInFieldSpeedRangePerTraversal;speedRangeInThisFieldLap];
                  
                  allLapInFieldSpikeTraversalAvgSpeeds=[allLapInFieldSpikeTraversalAvgSpeeds; inLapInFieldSpikeTraversalAvgSpeeds(:)];
                  allLapInFieldSpikeTraversalSpeedRanges=[allLapInFieldSpikeTraversalSpeedRanges; inLapInFieldSpikeTraversalSpeedRanges(:)];
                  
                  allLapInFieldSpikeTimes=[allLapInFieldSpikeTimes; inLapInFieldSpikeTimesInLap(:)];
                   
                   allLapInFieldSpikeTimesInExp=[allLapInFieldSpikeTimesInExp;inLapInFieldSpikeTimesInExp(:)];
                   allLapInFieldSpikeLapNums=[allLapInFieldSpikeLapNums; inLapInFieldSpikeLapNums(:)];
                   allLapInFieldSpikePhases=[allLapInFieldSpikePhases; inLapInFieldSpikePhases(:)];
                   allLapSpikeCycleIDsFromFieldEntry=[allLapSpikeCycleIDsFromFieldEntry; spikeCycleIDsFromFieldEntry(:)];
                   allLapInFieldCycleMinTimes=[allLapInFieldCycleMinTimes;inLapInFieldCycleMinTimes(:)];
                   allLapInFieldSpikeSpeeds=[allLapInFieldSpikeSpeeds;inLapInFieldSpikeSpeeds(:)];
                   
                   allLapInFieldSpikesInLapWithStop=[allLapInFieldSpikesInLapWithStop;inLapInFieldSpikesInLapWithStop(:)];
                   allLapInFieldSpikesInLapWithLateStart=[allLapInFieldSpikesInLapWithLateStart;inLapInFieldSpikesInLapWithLateStart(:)];
                   
                   allLapInFieldSpikeSpeedsSigned=[allLapInFieldSpikeSpeedsSigned;inLapInFieldSpikeSpeedsSigned(:)];
                   
		allLapInFieldSpikeDists=[allLapInFieldSpikeDists; inLapInFieldSpikeDists(:)];
		allLapInFieldSpikeDistsFieldFrac=[allLapInFieldSpikeDistsFieldFrac; inLapInFieldSpikeDists(:)/currFieldWidthM];
        
        allLapInFieldNumCycles=[allLapInFieldNumCycles;numCyclesInLapInField];
        
        allLapInFieldMeanCycleDurs=[allLapInFieldMeanCycleDurs;meanInFieldCycleDuration];

                %{
                 subplot(2,10,20)
                  plot(inLapInFieldSpikeDists/fieldWidth,inLapInFieldSpikePhases,'k.','MarkerSize',10)
                    hold on
                 %if(li==1)
                hold on
                xlim([0 1])
                 ylim([0 360])
             
               
                   xlabel('Dist (frac field)')
                 yticklabels({})
                 %ylabel('Theta phase (deg)')
                 title('Dist vs phase')
                %}
                 %end
                
                  
                
                %currLapRates=gaussSpikeRatePerTime(inLapIdxes);
                %durLength=min(length(currLapRates),numTimePts);
             
                %gaussRateProfile(li,1:durLength)=currLapRates(1:durLength);
                %set(gcf,'Visible','off')
                %}
                
            end %lap loop  
            
                if(showRaster)
                    if(isvalid(h1))
                    axes(h1)

                       %xlabel('Time since field entry (frac of avg field duration)')

                    ylabel('Lap number')
                    title('Temporal spike raster per lap')
                    ylim([-Inf Inf])
                    
                    end
               

                 figure(fRaster)
                uberTitle(removeUnderscores(sprintf('%s spike raster since field entry per lap',fileBaseName)))
         setFigFontTo(18)
                maxFig
                %plot(timeSinceAxis,nanmean(gaussRateProfile,1))
                saveas(gcf,sprintf('sinceFieldEntrySpikeRaster%s.tif',fileBaseName))
                end
                
           if(showPlots)
                figure(fPhaseVsCycleNum)
                ylim([0 360])
                uberTitle(removeUnderscores(sprintf('%s spike phase since field entry per lap',fileBaseName)))
                

         setFigFontTo(18)
                maxFig
                %plot(timeSinceAxis,nanmean(gaussRateProfile,1))
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %SAVE TIME VS PHASE IMAGE
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                saveas(gcf,saveImageFilePath)

                %{
                disp('')
                drawnow
                
                %}
                close all
           end
        
         thetaBasedTimeVsPhaseInfo.(currFieldDirection).(currFieldStr).allLapInFieldMeanPhasePerTraversal=allLapInFieldMeanPhasePerTraversal;
         thetaBasedTimeVsPhaseInfo.(currFieldDirection).(currFieldStr).allLapInFieldMeanSpeedPerTraversal=allLapInFieldMeanSpeedPerTraversal;
         thetaBasedTimeVsPhaseInfo.(currFieldDirection).(currFieldStr).allLapInFieldSpeedRangePerTraversal=allLapInFieldSpeedRangePerTraversal;
         
         thetaBasedTimeVsPhaseInfo.(currFieldDirection).(currFieldStr).allLapInFieldAllTimeVsPhasePerTraversal=allLapInFieldAllTimeVsPhasePerTraversal;
         
        thetaBasedTimeVsPhaseInfo.(currFieldDirection).(currFieldStr).numThetaCyclesInFieldPerLap=numThetaCyclesInFieldPerLap;
        thetaBasedTimeVsPhaseInfo.(currFieldDirection).(currFieldStr).allLapInFieldSpikeTimes=allLapInFieldSpikeTimes;
        thetaBasedTimeVsPhaseInfo.(currFieldDirection).(currFieldStr).allLapInFieldSpikeTimesInExp=allLapInFieldSpikeTimesInExp;
           thetaBasedTimeVsPhaseInfo.(currFieldDirection).(currFieldStr).allLapSpikeCycleIDsFromFieldEntry=allLapSpikeCycleIDsFromFieldEntry;
           thetaBasedTimeVsPhaseInfo.(currFieldDirection).(currFieldStr).allLapInFieldCycleMinTimes=allLapInFieldCycleMinTimes;
           thetaBasedTimeVsPhaseInfo.(currFieldDirection).(currFieldStr).allLapInFieldSpikeLapNums=allLapInFieldSpikeLapNums;
           
           thetaBasedTimeVsPhaseInfo.(currFieldDirection).(currFieldStr).allLapInFieldSpikeTraversalAvgSpeeds=allLapInFieldSpikeTraversalAvgSpeeds;
           %thetaBasedTimeVsPhaseInfo.(currFieldDirection).(currFieldStr).allLapInFieldSpikeTraversalSpeedRanges=allLapInFieldSpikeTraversalSpeedRanges;
           thetaBasedTimeVsPhaseInfo.(currFieldDirection).(currFieldStr).allLapInFieldSpikeTraversalSpeedStdDevs=allLapInFieldSpikeTraversalSpeedRanges;

           
           thetaBasedTimeVsPhaseInfo.(currFieldDirection).(currFieldStr).allLapInFieldSpikeSpeeds=allLapInFieldSpikeSpeeds;
           thetaBasedTimeVsPhaseInfo.(currFieldDirection).(currFieldStr).allLapInFieldSpikeSpeedsSigned=allLapInFieldSpikeSpeedsSigned;
           thetaBasedTimeVsPhaseInfo.(currFieldDirection).(currFieldStr).allLapInFieldSpikesInLapWithStop=allLapInFieldSpikesInLapWithStop;
           thetaBasedTimeVsPhaseInfo.(currFieldDirection).(currFieldStr).allLapInFieldSpikesInLapWithLateStart=allLapInFieldSpikesInLapWithLateStart;
           
               thetaBasedTimeVsPhaseInfo.(currFieldDirection).(currFieldStr).allLapInFieldSpikePhases=allLapInFieldSpikePhases;
           thetaBasedTimeVsPhaseInfo.(currFieldDirection).(currFieldStr).allLapInFieldSpikeDists=allLapInFieldSpikeDists;
           thetaBasedTimeVsPhaseInfo.(currFieldDirection).(currFieldStr).allLapInFieldSpikeDistsFieldFrac=allLapInFieldSpikeDistsFieldFrac;
            thetaBasedTimeVsPhaseInfo.(currFieldDirection).(currFieldStr).allLapInFieldNumCycles=allLapInFieldNumCycles;
            thetaBasedTimeVsPhaseInfo.(currFieldDirection).(currFieldStr).allLapInFieldMeanCycleDurs=allLapInFieldMeanCycleDurs;
            
            allSpikeGoodLapIdxes=~(allLapInFieldSpikesInLapWithStop & allLapInFieldSpikesInLapWithLateStart);
            
            thetaBasedTimeVsPhaseInfo.(currFieldDirection).(currFieldStr).allSpikeGoodLapIdxes=allSpikeGoodLapIdxes;
            

               numThetaCyclesInFieldPerLap=[];
                allLapInFieldSpikeTimes=[];
                allLapInFieldSpikeTimesInExp=[];
                allLapInFieldSpikeTraversalAvgSpeeds=[];
               allLapInFieldSpikePhases=[];
               allLapSpikeCycleIDsFromFieldEntry=[];
               allLapInFieldCycleMinTimes=[];
               allLapInFieldSpikeSpeeds=[];
               allLapInFieldSpikeSpeedsSigned=[];
               allLapInFieldSpikesInLapWithStop=[];
               allLapInFieldSpikesInLapWithLateStart=[];
               allLapInFieldSpikeDists=[];
               allLapInFieldSpikeDistsFieldFrac=[];
               allLapInFieldNumCycles=[];
               allLapInFieldMeanCycleDurs=[];
               
               allLapInFieldMeanPhasePerTraversal=[];
               allLapInFieldMeanSpeedPerTraversal=[];
               allLapInFieldSpeedRangePerTraversal=[];
               
            %close(fRaster)
            %close(fPhaseVsCycleNum)
        end
        
        %subplot(2,1,1)
        
         %subplot(2,10,[1 9])
     
        %subplot(2,1,2)
           %subplot(2,10,[11 19])
           %{
           if(isvalid(h2))
               axes(h2)

               xlabel('Distance since field entry (frac of avg field width)')
        
               
            ylabel('Lap number')
            title('Spatial spike raster per lap')
             ylim([-Inf Inf])
           end
           %}
        
       
        
    end
    
    thetaBasedTimeVsPhaseInfo.minAbsSpeedInLap=minAbsSpeedInLap;
    thetaBasedTimeVsPhaseInfo.maxFirstCycleInField=maxFirstCycleInField;
    thetaBasedTimeVsPhaseInfo.minLapNumSpikes=minLapNumSpikes;
    thetaBasedTimeVsPhaseInfo.eliminateBackwards=eliminateBackwards;
    thetaBasedTimeVsPhaseInfo.firstSpikeAsStartOfTimeField=firstSpikeAsStartOfTimeField;
    thetaBasedTimeVsPhaseInfo.minNumSpikesPerTraversal=minNumSpikesPerTraversal;
    
    save(currFilePath,'thetaBasedTimeVsPhaseInfo','-append')
    thetaBasedTimeVsPhaseInfo=[];
end

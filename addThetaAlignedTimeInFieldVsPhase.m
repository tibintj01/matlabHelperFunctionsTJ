close all; clear all; clc
dataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles';
%xxxxxsetTightSubplots_SpaceTime

filePaths=getFilePathsRegex(dataDir,'*mat');

colorCodeForSpeed=0;
colorCodeForPhase=1;
showPlots=0;

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
    
    [currSessionNum,currSessName]=getSessionNum(fileBaseName);
    
    if(isnan(currSessName))
        continue
    end
    
    
    if(~contains(currSessName,'atsby'))
        %continue
    end
    
     [leftwardRefUnitData,rightwardRefUnitData]=getRefUnitData(currSessName);
    
    data=load(currFilePath);
    allLapInFieldSpikeTimes=[];
               allLapInFieldSpikePhases=[];
               allLapSpikeCycleIDsFromFieldEntry=[];
               allLapInFieldCycleMinTimes=[];
               allLapInFieldSpikeSpeeds=[];
               allLapInFieldMeanCycleDurs=[];
               

    for di=1:2
            %figure(1)
            
            allLapInFieldSpikeTimes=[];
               allLapInFieldSpikePhases=[];
               allLapSpikeCycleIDsFromFieldEntry=[];
               allLapInFieldCycleMinTimes=[];
               allLapInFieldSpikeSpeeds=[];
        if(di==2)
            currFieldDirection='leftward'; 
            goingBackwardIdxes=data.speedPerTimeStep>0;
            gaussSpikeRatePerTime=data.leftLocalGaussSpikeRatePerPosTime;
            refUnitData=leftwardRefUnitData;
            
            if(useRefThetaCh && ~isnan(data.refThetaCh))
             spikePhases=data.leftSpikesRefPhasesMaxToMax;
            else
             spikePhases=data.leftSpikePhases;
   
            end
            
            spikeTimes=data.leftSpikeTimes;
            spikePositions=data.leftSpikePositions;

            if(justFirstSpikes)
                if(~isfield(data,'leftIsFirstSpikeInRefCycle'))
                    continue
                end
                spikeTimes=spikeTimes(data.leftIsFirstSpikeInRefCycle);
                spikePositions=spikePositions(data.leftIsFirstSpikeInRefCycle);
                spikePhases=spikePhases(data.leftIsFirstSpikeInRefCycle);
            end
            
        elseif(di==1)
            if(strContainsCircSessionName(currFileName))
                continue %circular mazes always leftward laps
            end
            currFieldDirection='rightward'; 
            goingBackwardIdxes=data.speedPerTimeStep<0;
            refUnitData=rightwardRefUnitData;
                      
            gaussSpikeRatePerTime=data.rightLocalGaussSpikeRatePerPosTime;
            
            if(useRefThetaCh && ~isnan(data.refThetaCh))
              spikePhases=data.rightSpikesRefPhasesMaxToMax;
            else
             spikePhases=data.rightSpikePhases;
         
            end
           
            spikeTimes=data.rightSpikeTimes;
            spikePositions=data.rightSpikePositions;

            if(justFirstSpikes)
                if(~isfield(data,'rightIsFirstSpikeInRefCycle'))
                    continue
                end
                spikeTimes=spikeTimes(data.rightIsFirstSpikeInRefCycle);
                spikePositions=spikePositions(data.rightIsFirstSpikeInRefCycle);
                spikePhases=spikePhases(data.rightIsFirstSpikeInRefCycle);
            end
        end
       
    spikePhasesPreShift=spikePhases;
     if(~useRefThetaCh || isnan(data.refThetaCh))
        spikePhases=manualApproxShiftDeg(spikePhasesPreShift,currSessName);
     else 
        spikePhases=manualApproxShiftRefDeg(spikePhasesPreShift,currSessName);
     end 
     
        
        spikeSpeeds=abs(interp1(data.positionTimeAxis,data.filledSignedSpeedPerTime,spikeTimes));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %load space time coordinates in field per time
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if(~isfield(data.directionSpecificStats,currFieldDirection))
            close all
            continue
        end

        numFields=size(data.directionSpecificStats.(currFieldDirection).expectedWidthPerField,1);
        
        lapStartTimes=refUnitData.lapStartTimesPerDir.(currFieldDirection);
        lapEndTimes=refUnitData.lapStopTimesPerDir.(currFieldDirection);
        
        refCycleInfo=load(refUnitData.unitInfo.cycleInfoFileName);
        refCycleMaxTimes=refCycleInfo.cycleMaxTimes;
        
        numLaps=length(lapStartTimes);
        
        positionTimeAxis=data.positionTimeAxis;
        
        refCycleMaxTimes=refCycleMaxTimes(refCycleMaxTimes>=positionTimeAxis(1) & refCycleMaxTimes<=positionTimeAxis(end));
        
        approxTimeStep=median(diff(positionTimeAxis));
        maxTimeSince=5;%seconds
        maxTimeSince=1.5;%seconds
         maxTimeSince=0.75;%seconds
         maxDistSince=1; %cm
         
        timeSinceAxis=(approxTimeStep/2):approxTimeStep:maxTimeSince;
        numTimePts=length(timeSinceAxis);
        
        gaussRateProfile=NaN(numLaps,numTimePts);
          
        currFieldStartPosPerField=data.(sprintf('manualFieldStartsM%s',currFieldDirection));
        currFieldEndPosPerField=data.(sprintf('manualFieldEndsM%s',currFieldDirection));

        for fii=1:numFields
            
            allLapInFieldSpikeTimes=[];
               allLapInFieldSpikePhases=[];
               allLapSpikeCycleIDsFromFieldEntry=[];
               allLapInFieldCycleMinTimes=[];
               allLapInFieldSpikeSpeeds=[];
               allLapInFieldSpikeDists=[];
               allLapInFieldSpikeDistsFieldFrac=[];
               allLapInFieldNumCycles=[];
               
               currFieldStr=sprintf('field%d',fii);
               
               if(showPlots)
                    fRaster=figure(20);
                    fPhaseVsCycleNum=figure(21);
                    autoArrangeFigures()
               end
            if(fii>length(data.directionSpecificStats.(currFieldDirection).avgFieldDurations))
                continue
            end
            
            try
                fieldDurLim=1.5*data.directionSpecificStats.(currFieldDirection).avgFieldDurations(fii);
                fieldEntryTimePerTime=data.directionSpecificStats.(currFieldDirection).fieldEntryTimePerTime(fii,:);
                fieldExitTimePerTime=data.directionSpecificStats.(currFieldDirection).fieldExitTimePerTime(fii,:);
                fieldWidth=abs(data.directionSpecificStats.(currFieldDirection).avgFieldWidths(fii));
                fieldDur=abs(data.directionSpecificStats.(currFieldDirection).avgFieldDurations(fii));
            catch
                continue
            end
            
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
            for li=1:numLaps
                inLapIdxes=positionTimeAxis>=lapStartTimes(li) & positionTimeAxis < lapEndTimes(li);
                fieldEntryTimeThisLap=nanmean(fieldEntryTimePerTime(inLapIdxes));
                fieldExitTimeThisLap=nanmean(fieldExitTimePerTime(inLapIdxes));
                
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
                
                inLapInFieldSpikeTimes=spikeTimes(inLapInFieldSpikeTimeIdxes);
                
               
                inLapInFieldSpikeSpeeds=spikeSpeeds(inLapInFieldSpikeTimeIdxes);
                inLapInFieldSpikePositions=spikePositions(inLapInFieldSpikeTimeIdxes);
                inLapInFieldSpikePhases=spikePhases(inLapInFieldSpikeTimeIdxes);
                
                %spikeCycleIDsFromFieldEntry=NaN(size(inLapInFieldSpikeTimes));
                
                if(length(inLapInFieldCycleMinTimes)>1)
                    spikeCycleIDsFromFieldEntry=floor(interp1(inLapInFieldCycleMinTimes,1:(length(inLapInFieldCycleMinTimes)),inLapInFieldSpikeTimes));
                else
                    spikeCycleIDsFromFieldEntry=1;
                end
               
                %inLapInFieldSpikeTimes=inLapInFieldSpikeTimes-min(inLapInFieldSpikeTimes);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %zero relative to theta cycle timing so that theta is time
                %keeper across laps
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                inLapInFieldSpikeTimes=inLapInFieldSpikeTimes-closestThetaCycleMinTime;
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
                inLapInFieldSpikeTimes(badSpeedIdxes)=[];
                inLapInFieldSpikeDists(badSpeedIdxes)=[];
                inLapInFieldSpikePhases(badSpeedIdxes)=[];
                
                try
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
                end

                if(colorCodeForSpeed)
                    speedIdxes=ceil(inLapInFieldSpikeSpeeds*100);
                    speedIdxes(speedIdxes>100)=100;
                    speedIdxes(isnan(speedIdxes))=100;
                    tickColor=speedColors(speedIdxes,:);
                elseif(colorCodeForPhase)
                    phaseIdxes=ceil(inLapInFieldSpikePhases);
                    phaseIdxes(phaseIdxes==0)=1; %ceil for exact 0
                    tickColor=phaseColors(phaseIdxes,:);
                end
                
                 if(isempty(inLapInFieldSpikeTimes))
                    continue
                end
                %subplot(2,1,1)
                if(showPlots)
                    figure(fPhaseVsCycleNum)
                    plot(spikeCycleIDsFromFieldEntry,inLapInFieldSpikePhases,'k.','MarkerSize',10)
                    hold on
                   xlabel('Theta cycle since field entry (num)')
                      ylabel('Theta phase (deg)')

                      xlim([0 maxTimeDisp/0.1])


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
                
                %plotRasterStyle(inLapInFieldSpikeTimes/fieldDur,li,NaN,NaN,tickColor);
                inLapInFieldSpikePhases=inLapInFieldSpikePhases(inLapInFieldSpikeTimes<=maxTimeDisp);
                inLapInFieldSpikeSpeeds=inLapInFieldSpikeSpeeds(inLapInFieldSpikeTimes<=maxTimeDisp);
                spikeCycleIDsFromFieldEntry=spikeCycleIDsFromFieldEntry(inLapInFieldSpikeTimes<=maxTimeDisp);
                inLapInFieldSpikeDists=inLapInFieldSpikeDists(inLapInFieldSpikeTimes<=maxTimeDisp);
               
		inLapInFieldSpikeTimes=inLapInFieldSpikeTimes(inLapInFieldSpikeTimes<=maxTimeDisp);
		 
                inLapInFieldCycleMinTimes=inLapInFieldCycleMinTimes(inLapInFieldCycleMinTimes<=maxTimeDisp);
                numCyclesInLapInField=length(inLapInFieldCycleMinTimes);
                
                if(showPlots)
                plotRasterStyle(inLapInFieldSpikeTimes,li,NaN,NaN,tickColor);
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
                  
                  if(showPlots)
                      
                      subplot(10,1,[7 8 9 10])
                    %plot(inLapInFieldSpikeTimes/fieldDur,inLapInFieldSpikePhases,'k.','MarkerSize',10)

                   plot(inLapInFieldSpikeTimes,inLapInFieldSpikePhases,'k.','MarkerSize',10)


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
                  
                  
                   allLapInFieldSpikeTimes=[allLapInFieldSpikeTimes; inLapInFieldSpikeTimes(:)];
                   allLapInFieldSpikePhases=[allLapInFieldSpikePhases; inLapInFieldSpikePhases(:)];
                   allLapSpikeCycleIDsFromFieldEntry=[allLapSpikeCycleIDsFromFieldEntry; spikeCycleIDsFromFieldEntry(:)];
                   allLapInFieldCycleMinTimes=[allLapInFieldCycleMinTimes;inLapInFieldCycleMinTimes(:)];
                   allLapInFieldSpikeSpeeds=[allLapInFieldSpikeSpeeds;inLapInFieldSpikeSpeeds(:)];
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
                
                  
                
                currLapRates=gaussSpikeRatePerTime(inLapIdxes);
                durLength=min(length(currLapRates),numTimePts);
             
                gaussRateProfile(li,1:durLength)=currLapRates(1:durLength);
                %set(gcf,'Visible','off')
                %}
                
            end %lap loop  
            
                if(showPlots)
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

                figure(fPhaseVsCycleNum)
                uberTitle(removeUnderscores(sprintf('%s spike phase since field entry per lap',fileBaseName)))
                ylim([0 360])

         setFigFontTo(18)
                maxFig
                %plot(timeSinceAxis,nanmean(gaussRateProfile,1))
                saveas(gcf,sprintf('sinceFieldEntrySpikePhaseVsTheta%s.tif',fileBaseName))

                disp('')
                drawnow
                close all
         end
        
        thetaBasedTimeVsPhaseInfo.(currFieldDirection).(currFieldStr).numThetaCyclesInFieldPerLap=numThetaCyclesInFieldPerLap;
        thetaBasedTimeVsPhaseInfo.(currFieldDirection).(currFieldStr).allLapInFieldSpikeTimes=allLapInFieldSpikeTimes;
           thetaBasedTimeVsPhaseInfo.(currFieldDirection).(currFieldStr).allLapSpikeCycleIDsFromFieldEntry=allLapSpikeCycleIDsFromFieldEntry;
           thetaBasedTimeVsPhaseInfo.(currFieldDirection).(currFieldStr).allLapInFieldCycleMinTimes=allLapInFieldCycleMinTimes;
           
           
           thetaBasedTimeVsPhaseInfo.(currFieldDirection).(currFieldStr).allLapInFieldSpikeSpeeds=allLapInFieldSpikeSpeeds;
               thetaBasedTimeVsPhaseInfo.(currFieldDirection).(currFieldStr).allLapInFieldSpikePhases=allLapInFieldSpikePhases;
           thetaBasedTimeVsPhaseInfo.(currFieldDirection).(currFieldStr).allLapInFieldSpikeDists=allLapInFieldSpikeDists;
           thetaBasedTimeVsPhaseInfo.(currFieldDirection).(currFieldStr).allLapInFieldSpikeDistsFieldFrac=allLapInFieldSpikeDistsFieldFrac;
            thetaBasedTimeVsPhaseInfo.(currFieldDirection).(currFieldStr).allLapInFieldNumCycles=allLapInFieldNumCycles;
            thetaBasedTimeVsPhaseInfo.(currFieldDirection).(currFieldStr).allLapInFieldMeanCycleDurs=allLapInFieldMeanCycleDurs;

               numThetaCyclesInFieldPerLap=[];
                allLapInFieldSpikeTimes=[];
               allLapInFieldSpikePhases=[];
               allLapSpikeCycleIDsFromFieldEntry=[];
               allLapInFieldCycleMinTimes=[];
               allLapInFieldSpikeSpeeds=[];
               allLapInFieldSpikeDists=[];
               allLapInFieldSpikeDistsFieldFrac=[];
               allLapInFieldNumCycles=[];
               allLapInFieldMeanCycleDurs=[];
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
    
    save(currFilePath,'thetaBasedTimeVsPhaseInfo','-append')
end

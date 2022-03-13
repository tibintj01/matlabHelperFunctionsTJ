close all; clear all; clc
dataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles';
%xxxxxsetTightSubplots_SpaceTime

filePaths=getFilePathsRegex(dataDir,'*mat');

colorCodeForSpeed=0;
colorCodeForPhase=1;

justFirstSpikes=1;
justFirstSpikes=0;

useRefThetaCh=0;
useRefThetaCh=1;


%autoShiftOffset=1;

cmapName=magma(100);
%cmapName=plasma(100);
speedColors=cmapName;

phaseColors=hsv(360);
%phaseColors=phasemap(360);

for fi=1:length(filePaths)
    %for fi=76:length(filePaths)
       fi
    currFilePath=filePaths{fi};
    currFileName=getFileNameFromPath(currFilePath);
    fileBaseName=currFileName(1:(end-4));
    
    [currSessionNum,currSessName]=getSessionNum(fileBaseName);
    
    if(isnan(currSessName))
        continue
    end
    
    %{
    if(~contains(currSessName,'uddy'))
        continue
    end
    %}
     [leftwardRefUnitData,rightwardRefUnitData]=getRefUnitData(currSessName);
    
    data=load(currFilePath);

    for di=1:2
            %figure(1)
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
                spikeTimes=spikeTimes(data.rightIsFirstSpikeInRefCycle);
                spikePositions=spikePositions(data.rightIsFirstSpikeInRefCycle);
                spikePhases=spikePhases(data.rightIsFirstSpikeInRefCycle);
            end
        end
       
    spikePhasesPreShift=spikePhases;
     if(~useRefThetaCh || isnan(data.refThetaCh))
        spikePhases=manualApproxShiftDeg(spikePhasesPreShift,currSessName);
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
        numLaps=length(lapStartTimes);
        
        positionTimeAxis=data.positionTimeAxis;
        
        approxTimeStep=median(diff(positionTimeAxis));
        maxTimeSince=5;%seconds
        maxTimeSince=1.5;%seconds
         maxTimeSince=0.75;%seconds
         maxDistSince=1; %cm
         
        timeSinceAxis=(approxTimeStep/2):approxTimeStep:maxTimeSince;
        numTimePts=length(timeSinceAxis);
        
        gaussRateProfile=NaN(numLaps,numTimePts);
           figure;
        currFieldStartPosPerField=data.(sprintf('manualFieldStartsM%s',currFieldDirection));

        for fii=1:numFields
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
            
            for li=1:numLaps
                inLapIdxes=positionTimeAxis>=lapStartTimes(li) & positionTimeAxis < lapEndTimes(li);
                fieldEntryTimeThisLap=nanmean(fieldEntryTimePerTime(inLapIdxes));
                fieldExitTimeThisLap=nanmean(fieldExitTimePerTime(inLapIdxes));
                
                inLapSpikeTimeIdxes=spikeTimes>=lapStartTimes(li) & spikeTimes < lapEndTimes(li);
                inLapInFieldSpikeTimeIdxes= inLapSpikeTimeIdxes & spikeTimes>=fieldEntryTimeThisLap & spikeTimes < fieldExitTimeThisLap;
               
                
                if(max(inLapInFieldSpikeTimeIdxes)==0)
                    continue
                end
                inLapInFieldSpikeTimes=spikeTimes(inLapInFieldSpikeTimeIdxes);
                inLapInFieldSpikeSpeeds=spikeSpeeds(inLapInFieldSpikeTimeIdxes);
                inLapInFieldSpikePositions=spikePositions(inLapInFieldSpikeTimeIdxes);
                inLapInFieldSpikePhases=spikePhases(inLapInFieldSpikeTimeIdxes);
                
                inLapInFieldSpikeTimes=inLapInFieldSpikeTimes-min(inLapInFieldSpikeTimes);
                %inLapInFieldSpikeDists=abs(inLapInFieldSpikePositions-min(inLapInFieldSpikePositions));
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
                
                
                %subplot(2,1,1)
     
               
                if(li==1 || ~isvalid(h1))
                 h1=subplot(2,10,1:9);
                else
                    hold(h1,'on')
                     axes(h1)
                end
                plotRasterStyle(inLapInFieldSpikeTimes/fieldDur,li,NaN,NaN,tickColor);
                colormap(gca,magma)
               
                
                cb=colorbar;
                ylabel(cb,'running speed at spike time (cm/s)')
               
                %xlim([0 maxTimeSince])
                xlim([0 1])
                caxis([0 100])
                 if(colorCodeForPhase)
                     colormap(gca,hsv)
                     %colormap(gca,phasemap)
                     if(useRefThetaCh)
                           ylabel(cb,'Ref. theta phase (deg)')
                       else
                         ylabel(cb,'Local theta phase (deg)')
                       end
                    caxis([0 360])
                end
                %drawnow
                
              
                
                 %subplot(2,1,2)
                
              
                 %subplot(2,10,11:19)
                 if(li==1 || ~isvalid(h2))
                    h2=subplot(2,10,11:19);
                else
                    hold(h2,'on')
                     axes(h2)
                end
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
                
                
                
                  subplot(2,10,10)
                plot(inLapInFieldSpikeTimes/fieldDur,inLapInFieldSpikePhases,'k.','MarkerSize',10)
                hold on
                %uberTitle(removeUnderscores(currSessName))
                %if(li==1)
                xlim([0 1])
                 ylim([0 360])
               
                    xlabel('Time (frac field)')
                  %ylabel('Theta phase (deg)')
                    yticklabels({})
                    title('Time vs phase')
                %end
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
                 %end
                
                  
                currLapRates=gaussSpikeRatePerTime(inLapIdxes);
                durLength=min(length(currLapRates),numTimePts);
             
                gaussRateProfile(li,1:durLength)=currLapRates(1:durLength);
                set(gcf,'Visible','off')

            end  
        end
        
        %subplot(2,1,1)
        
         %subplot(2,10,[1 9])
         if(isvalid(h1))
             axes(h1)

               xlabel('Time since field entry (frac of avg field duration)')
          
            ylabel('Lap number')
            title('Temporal spike raster per lap')
            ylim([-Inf Inf])
         end
        %subplot(2,1,2)
           %subplot(2,10,[11 19])
           if(isvalid(h2))
               axes(h2)

               xlabel('Distance since field entry (frac of avg field width)')
        
               
            ylabel('Lap number')
            title('Spatial spike raster per lap')
             ylim([-Inf Inf])
           end
        
       
        uberTitle(removeUnderscores(sprintf('%s spike raster since field entry per lap',fileBaseName)))
 setFigFontTo(14)
        maxFig
        %plot(timeSinceAxis,nanmean(gaussRateProfile,1))
        saveas(gcf,sprintf('sinceFieldEntrySpikeRaster%s.tif',fileBaseName))
        disp('')
        close all
        
    end
    
end
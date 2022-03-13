close all; clear all; clc
dataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles';
cycleDataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/';
showPlots=1;
%setTightSubplots_SpaceTime

%Cicero_09172014 has only 1 strongly precessing cell used
%sessionNames={'Achilles_10252013','Achilles_11012013','Buddy_06272013', 'Cicero_09102014','Cicero_09172014','Gatsby_08282013'};
sessionNames={'Achilles_10252013','Achilles_11012013','Buddy_06272013', 'Cicero_09102014','Gatsby_08282013'};

sessionNames={'Achilles_11012013'};
%chosenUnitNums=[1024 1022 107 407];
%chosenFieldIDs=[1 1 3 2];

chosenUnitNums=[407 1024 1022 409 308];
chosenUnitNums=[ 1024 1022 409 308];
chosenFieldIDs=[1 1 1 2];
onlyChosen=1;

%sessionNames={'Gatsby_08282013'};

halfLineHeight=0.01; %raster line height

lowSpeedMax=0.3;
highSpeedMin=0.6;

useReferenceFieldEntryTime=0;
useTimeInLap=0;
%useTimeInLap=1;
twoSpeeds=1;
twoSpeeds=0;

timeMax=2.5;
timeMax=2;
%timeMax=3;
timeMax=2;

%{
lowSpeedMax=0.35;
highSpeedMin=0.7;

lowSpeedMax=0.25;
highSpeedMin=0.5;
%}
lowSpeedMax=0.4;
highSpeedMin=0.5;

lowSpeedMin=0.05;
lowSpeedMax=0.3;
highSpeedMin=0.6;

lowSpeedMin=0.05;
highSpeedMax=1.2;
lowSpeedMax=0.35;
highSpeedMin=0.7;
lowSpeedMax=0.3;
highSpeedMin=0.6;
%{
lowSpeedMax=0.4;
highSpeedMin=0.8;
%}


%{
lowSpeedMax=0.1;
highSpeedMin=1;
%}

%{
lowSpeedMax=0.5;
highSpeedMin=0.5;
%}

%{
lowSpeedMax=1.5;
highSpeedMin=0;
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%loop through each session and all corresponding units
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for si=1:length(sessionNames) 
    currSessName=sessionNames{si};
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %separate lap nums into top and bottom 20 percentile average speed
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      [minTime, maxTime,numLapsPerDir]=getSessionTimeBounds(sessionNames{si});
    %maxTime=minTime+10;
    
    [leftwardRefUnitData,rightwardRefUnitData]=getRefUnitData(currSessName);
    
    currSessUnitDataPaths=getRegexFilePaths(dataDir,sprintf('%s*Unit*.mat',currSessName));
        numUnitsInSess=length(currSessUnitDataPaths);
        
        for di=1:2
                tic
                if(di==1)
                    currDir='right';
                    currDirStr='rightward';
                else
                    currDir='left';
                    currDirStr='leftward';
                end

                for fi=1:numUnitsInSess
                    %currUnitData=currSessMasterUnitData{fi};
                    currUnitData=load(currSessUnitDataPaths{fi});
                    unitIDnum=currUnitData.unitInfo.unitIDnum;
                    
                    if(onlyChosen)
                        if(isempty(intersect(unitIDnum,chosenUnitNums)))
                            continue
                        end
                    end
                   
                    spikeTimes=currUnitData.unitInfo.spikeTimes;
                    positionTimeAxis=currUnitData.positionTimeAxis;
                    positionPerTimeStep=currUnitData.positionPerTimeStep;
                    
                    unitMinPos=min(positionPerTimeStep);
                    unitMaxPos=max(positionPerTimeStep);
                    
                    timeInLapPerTimeStep=NaN(size(positionPerTimeStep));
                    
                    
                    if(~isfield(currUnitData.lapStartTimesPerDir,currDirStr))
                        continue
                    end
                    
                    
                    lapStartTimesPerDir=currUnitData.lapStartTimesPerDir.(currDirStr);
                    lapStopTimesPerDir=currUnitData.lapStopTimesPerDir.(currDirStr);
                    numLaps=length(lapStartTimesPerDir);
                    
                    entryTimePerDirPerFieldPerLap=currUnitData.entryTimePerDirPerFieldPerLap;
                    exitTimePerDirPerFieldPerLap=currUnitData.exitTimePerDirPerFieldPerLap;
                    
                    
                     if(di>size(entryTimePerDirPerFieldPerLap,1))
                        continue
                     end
                      
                    currUnitFieldStartsPos=currUnitData.(sprintf('manualFieldStartsM%s',currDirStr));
                    currUnitFieldEndsPos=currUnitData.(sprintf('manualFieldEndsM%s',currDirStr));
                    numFieldsInDir=length(currUnitFieldStartsPos);
                    

                    timeInFieldPerTimeStep=NaN(numFieldsInDir,length(positionPerTimeStep));
                    
                    for fieldIdx=1:numFieldsInDir
                        
                        if(onlyChosen)
                            if(unitIDnum==chosenUnitNums(1) && fieldIdx ~=chosenFieldIDs(1))
                                continue
                            end

                            if(unitIDnum==chosenUnitNums(2) && fieldIdx ~=chosenFieldIDs(2))
                                continue
                            end

                            if(unitIDnum==chosenUnitNums(3) && fieldIdx ~=chosenFieldIDs(3))
                                continue
                            end

                            if(unitIDnum==chosenUnitNums(4) && fieldIdx ~=chosenFieldIDs(4))
                                continue
                            end
                        end
                        
                        
                        currFieldPosStart=currUnitFieldStartsPos(fieldIdx);
                        currFieldPosEnd=currUnitFieldEndsPos(fieldIdx);
                        
                        currFieldPosHigh=max([currFieldPosStart currFieldPosEnd]);
                        currFieldPosLow=min([currFieldPosStart currFieldPosEnd]);

                        fieldWidth=currFieldPosHigh-currFieldPosLow;
                        currFieldPosHigh=currFieldPosHigh+fieldWidth*0.4;
                        currFieldPosLow=currFieldPosLow-fieldWidth*0.4;
    
                        
                        if(fieldIdx>size(entryTimePerDirPerFieldPerLap,2))
                            continue
                        end
                        
                        for ti=1:length(positionTimeAxis)
                            currTime=positionTimeAxis(ti);
                            [~,closestLapStartIdx]=min(abs(currTime-lapStartTimesPerDir));
                             %[~,closestFieldStartIdx]=min(abs(currTime-squeeze(entryTimePerDirPerFieldPerLap(di,fi,:))));

                            closestFieldStartIdx=closestLapStartIdx;
                            if(currTime<lapStartTimesPerDir(closestLapStartIdx))
                                closestLapStartIdx=closestLapStartIdx-1;
                            end
                            if(closestLapStartIdx<1 || closestLapStartIdx > numLaps)
                                continue
                            end
                              thisLapStartTime=lapStartTimesPerDir(closestLapStartIdx);
                            thisLapEndTime=lapStopTimesPerDir(closestLapStartIdx);
                             if(currTime<=thisLapEndTime)
                                timeInLapPerTimeStep(ti)=currTime-thisLapStartTime;
                             end

                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %field entry time calculations
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                            while(closestFieldStartIdx>size(entryTimePerDirPerFieldPerLap,3))
                                closestFieldStartIdx=closestFieldStartIdx-1;
                            end

                            if(currTime<entryTimePerDirPerFieldPerLap(di,fieldIdx,closestFieldStartIdx))
                                closestFieldStartIdx=closestFieldStartIdx-1;
                            end

                            if(closestFieldStartIdx<1 || closestFieldStartIdx > numLaps)
                                continue
                            end

                            %closestLapStartIdx=closestLapStartIdx-1;
                            %{
                             if(closestLapStartIdx<1 || closestLapStartIdx > numLaps)
                                continue
                            end
                            %}

                            thisLapFieldStartTime=entryTimePerDirPerFieldPerLap(di,fieldIdx,closestFieldStartIdx);
                            thisLapFieldEndTime=exitTimePerDirPerFieldPerLap(di,fieldIdx,closestFieldStartIdx);

                            if(currTime<=thisLapFieldEndTime)
                                timeInFieldPerTimeStep(fieldIdx,ti)=currTime-thisLapFieldStartTime;
                            end
                        end
                   
                    filledSignedSpeedPerTime=currUnitData.filledSignedSpeedPerTime;
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %filter spike times by speed at that time
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    spikeSpeeds=abs(interp1(positionTimeAxis,filledSignedSpeedPerTime,spikeTimes));
                    
                    %if(twoSpeeds)
                    %    lowSpeedSpikeIdxes=spikeSpeeds<=highSpeedMin & spikeSpeeds>=lowSpeedMin;
                    %else
                        lowSpeedSpikeIdxes=spikeSpeeds<=lowSpeedMax & spikeSpeeds>=lowSpeedMin;
                    %end
                    
                    highSpeedSpikeIdxes=spikeSpeeds<=highSpeedMax & spikeSpeeds>=highSpeedMin;
                    
                    if(twoSpeeds)
                        midSpeedSpikeIdxes=spikeSpeeds<-Inf;
                    else
                        midSpeedSpikeIdxes=spikeSpeeds>lowSpeedMax & spikeSpeeds<highSpeedMin;
                    end
                    
                    lowSpeedSpikeTimes=spikeTimes(lowSpeedSpikeIdxes);
                    highSpeedSpikeTimes=spikeTimes(highSpeedSpikeIdxes);
                    midSpeedSpikeTimes=spikeTimes(midSpeedSpikeIdxes);
                    
                    %{
                    [positionBins,lowSpeedFiringRatePerPositionRight, lowSpeedFiringRatePerPositionLeft]...
                        = getFiringRatePerPosition(lowSpeedSpikeTimes,positionTimeAxis,positionPerTimeStep,filledSignedSpeedPerTime);
                    
                    [positionBins,highSpeedFiringRatePerPositionRight, highSpeedFiringRatePerPositionLeft]...
                        = getFiringRatePerPosition(highSpeedSpikeTimes,positionTimeAxis,positionPerTimeStep,filledSignedSpeedPerTime);
                    
                    %}
                    
                  
                    timeInMetaFieldPerTimeStep=timeInFieldPerTimeStep(fieldIdx,:);
                      if(unitIDnum==1024)
                        save('refTimeInFieldPerTimeStep1024.mat','timeInMetaFieldPerTimeStep')
                      end
                    
                      if(useReferenceFieldEntryTime)
                          refData=load('refTimeInFieldPerTimeStep1024.mat');
                          timeInMetaFieldPerTimeStep=refData.timeInMetaFieldPerTimeStep;
                      end
                      if(useTimeInLap)
                         timeInMetaFieldPerTimeStep=timeInLapPerTimeStep;
                      end
                   
                    
                    [positionBins,lowSpeedFiringRatePerPositionRight, lowSpeedFiringRatePerPositionLeft,...
                        timeInLapBins,lowSpeedFiringRatePerTimeRight, lowSpeedFiringRatePerTimeLeft] ...
                        = getFiringRatePerPositionAndTime(lowSpeedSpikeTimes,positionTimeAxis,positionPerTimeStep,timeInMetaFieldPerTimeStep,filledSignedSpeedPerTime,currFieldPosLow,currFieldPosHigh);
                        %= getFiringRatePerPositionAndTime(lowSpeedSpikeTimes,positionTimeAxis,positionPerTimeStep,timeInLapPerTimeStep,filledSignedSpeedPerTime);

                    [positionBins,highSpeedFiringRatePerPositionRight, highSpeedFiringRatePerPositionLeft,...
                        timeInLapBins,highSpeedFiringRatePerTimeRight, highSpeedFiringRatePerTimeLeft] ...
                        = getFiringRatePerPositionAndTime(highSpeedSpikeTimes,positionTimeAxis,positionPerTimeStep,timeInMetaFieldPerTimeStep,filledSignedSpeedPerTime,currFieldPosLow,currFieldPosHigh);
                        %= getFiringRatePerPositionAndTime(highSpeedSpikeTimes,positionTimeAxis,positionPerTimeStep,timeInLapPerTimeStep,filledSignedSpeedPerTime);

                     [positionBins,midSpeedFiringRatePerPositionRight, midSpeedFiringRatePerPositionLeft,...
                        timeInLapBins,midSpeedFiringRatePerTimeRight, midSpeedFiringRatePerTimeLeft] ...
                        = getFiringRatePerPositionAndTime(midSpeedSpikeTimes,positionTimeAxis,positionPerTimeStep,timeInMetaFieldPerTimeStep,filledSignedSpeedPerTime,currFieldPosLow,currFieldPosHigh);

                    if(di==1)
                        subplot(2,1,1)
                        plot(positionBins,lowSpeedFiringRatePerPositionRight,'b-','LineWidth',4)
                        hold on
                        plot(positionBins,midSpeedFiringRatePerPositionRight,'k-','LineWidth',4)
                        plot(positionBins,highSpeedFiringRatePerPositionRight,'r-','LineWidth',4)
                        
                         subplot(2,1,2)
                        plot(timeInLapBins,lowSpeedFiringRatePerTimeRight,'b-','LineWidth',4)
                        hold on
                        plot(timeInLapBins,midSpeedFiringRatePerTimeRight,'k-','LineWidth',4)
                        plot(timeInLapBins,highSpeedFiringRatePerTimeRight,'r-','LineWidth',4)
                        
                    else
                     
                  
                        subplot(2,1,1)
                        plot(positionBins,lowSpeedFiringRatePerPositionLeft,'b-','LineWidth',4)
                        hold on
                        plot(positionBins,midSpeedFiringRatePerPositionLeft,'k-','LineWidth',4)
                        plot(positionBins,highSpeedFiringRatePerPositionLeft,'r-','LineWidth',4)
                     

                        subplot(2,1,2)
                        plot(timeInLapBins,lowSpeedFiringRatePerTimeLeft,'b-','LineWidth',4)
                        hold on
                          plot(timeInLapBins,midSpeedFiringRatePerTimeLeft,'k-','LineWidth',4)
                        plot(timeInLapBins,highSpeedFiringRatePerTimeLeft,'r-','LineWidth',4)
                      
                    end
                     
                     subplot(2,1,2)
                     if(useTimeInLap)
                        xlabel('Time in lap (sec)')
                     else
                        xlabel('Time in field (sec)')
                     end
                     ylabel('Firing rate (Hz)')
                     %legend(sprintf('%.2f<speed<=%.2f m/s',lowSpeedMin,lowSpeedMax),sprintf('%.2f m/s<speed<%.2f m/s',lowSpeedMax,highSpeedMin),sprintf('%.2f<=speed<%.2f m/s',highSpeedMin,highSpeedMax))
                     %axis tight
                     if(unitIDnum==308 || unitIDnum==407)
                         xlim([0 3])
                     else
                        xlim([0 timeMax])
                     end
                     title('Time warp and firing vs time')
                     
                     subplot(2,1,1)
                     xlabel('Position in track (m)')
                     ylabel('Firing rate (Hz)')
                     %legend(sprintf('%.2f<=speed<=%.2f m/s',lowSpeedMin,lowSpeedMax),sprintf('%.2f m/s<speed<%.2f m/s',lowSpeedMax,highSpeedMin),sprintf('speed>=%.2f m/s',highSpeedMin))
                      xlim([unitMinPos unitMaxPos])
                      title('Time warp and firing vs position')
                      
                      for si=1:2
                          subplot(2,1,si)
                          legend(sprintf('%.2f m/s<speed<%.2f m/s',lowSpeedMin,lowSpeedMax),sprintf('%.2f m/s<speed<%.2f m/s',lowSpeedMax,highSpeedMin),sprintf('%.2f m/s<speed<%.2f m/s',highSpeedMin,highSpeedMax))
                           box off
                      end
                      
                     uberTitle(removeUnderscores(sprintf('%s_unit_%d_Field%d_%s',currSessName,unitIDnum,fieldIdx,currDirStr)))
                     setFigFontTo(18)
                     maxFig
                    drawnow
                    saveas(gcf,sprintf('%s_unit_%d_Field%d_%s_HighVsLowSpeedFiring.tif',currSessName,unitIDnum, fieldIdx,currDirStr))
                    close all
                     end
                end
        end
end

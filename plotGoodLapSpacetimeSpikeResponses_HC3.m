
close all; clear all; clc
%dataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles';

processedDataDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData/';
unitDataDir=fullfile(processedDataDir,'unitSpikeInfoStructs');
dataDir=unitDataDir;

savePhaseJumpVsTimeDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3Images/phaseJumpVsTimeInField';
touchDir(savePhaseJumpVsTimeDir)

spaceTimeResponsesDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3Images/spaceTimeResponses';
touchDir(spaceTimeResponsesDir)

%GAUSSKERNELWIDTHSEC=0.1;
GAUSSKERNELWIDTHSEC=0.01;

maxSpeedRangeInField=0.1; %10 cm/s
%maxSpeedRangeInField=0.05; %10 cm/s
maxSpeedRangeInField=0.2; %20 cm/s
maxSpeedRangeInField=0.15; %15 cm/s

maxSpeedRangeInField=0.15; %15 cm/s
maxSpeedRangeInField=0.1; %10 cm/s

filePaths=getFilePathsRegex(dataDir,'*mat');
totalNumFields=0;

showPlots=0;
showPlots=1;

offsetAngle=180;
noCirc=0;
justAchilles=0;
zeroVal=0.5;
zeroVal=0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set shift!!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
autoShiftOffset=0;
autoShiftOffset=1;

if(~exist('masterSpacetimePhaseResponse.mat','file'))

masterDistPerTime=[];
masterLapTimePerTime=[];
masterSpikeRatePerTime=[];
masterSpikePhasePerTime=[];
masterSpeedPerTime=[];
masterFieldNumPerTime=[];
masterSesNumPerTime=[];
masterDirectionPerTime=[];


heatRes=0.005
heatRes=0.0075
heatRes=0.01

heatResUnit=0.02;
heatResUnit=0.05;

fieldCount=0;
cellCount=0;

startFileID=1;

totalUsedFieldCount=0;

for fi=startFileID:length(filePaths)
       fi
    currFilePath=filePaths{fi};
    

 
    currFileName=getFileNameFromPath(currFilePath);
    fileBaseName=currFileName(1:(end-4));
    
    %{
    [currSessionNum, currSessName]=getSessionNum(fileBaseName);
    [refUnitDataLeftward,refUnitDataRightward] = getRefUnitData(currSessName);
    checkSpeed(currFilePath);
    %}
    
    data=load(currFilePath);
    
    posInfo=load(data.unitInfo.positionInfoForSessionSavePath);
    speedPerTimeStep=posInfo.speedAlongTrackPerTimeMSec;
    
    spikeTimeDirectionAssignment=data.unitInfo.spikeTimeDirectionAssignment;
    
    

    for di=1:2
            %figure(1)
        if(di==2)
            currFieldDirection='leftward'; 
            goingBackwardIdxes=speedPerTimeStep>0;
            
            refUnitData=data;
            
            
             %unitDirSpikeTimes=data.leftSpikeTimes;
             unitDirSpikeTimes=data.unitInfo.spikeTimes(spikeTimeDirectionAssignment==2);
            unitDirSpikeRefPhases=data.unitInfo.leftSpikePhases;
            unitDirSpikeRefCycleIDs=data.leftSpikeCycleIDs;
          
            
        elseif(di==1)
            currFieldDirection='rightward'; 
            goingBackwardIdxes=speedPerTimeStep<0;
            
            refUnitData=data;

              %unitDirSpikeTimes=data.rightSpikeTimes;
              unitDirSpikeTimes=data.unitInfo.spikeTimes(spikeTimeDirectionAssignment==1);
              %{
              if(~isfield(data,'rightSpikesRefPhases'))
                  continue
              end
              %}
            unitDirSpikeRefPhases=data.unitInfo.rightSpikePhases;
            unitDirSpikeRefCycleIDs=data.rightSpikeCycleIDs;
        end
        
        if(isfield(refUnitData,'lapStartTimesPerDir'))
            
            if(isfield(refUnitData.lapStartTimesPerDir,currFieldDirection))
                lapStartTimes=refUnitData.lapStartTimesPerDir.(currFieldDirection);
                lapStopTimes=refUnitData.lapStopTimesPerDir.(currFieldDirection);
            else
                continue
            end
            
        else
            continue
        end
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %load space time coordinates in field per time
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if(~isfield(data.directionSpecificStats,currFieldDirection))
            close all
            continue
        end
        localDistInFieldPerTime=data.directionSpecificStats.(currFieldDirection).localDistInFieldPerTime;
        localTimeInFieldPerTime=data.directionSpecificStats.(currFieldDirection).localTimeInFieldPerTime;
        
        localNormDistInFieldPerTime=data.directionSpecificStats.(currFieldDirection).localNormDistInFieldPerTime;
        localNormTimeInFieldPerTime=data.directionSpecificStats.(currFieldDirection).localNormTimeInFieldPerTime;
        
        localSpikePhasePerTime=data.directionSpecificStats.(currFieldDirection).localSpikePhasePerTime;
        %if(isfield(data,'localGaussSpikeRatePerPosTime'))
        if(di==1)
            localGaussSpikeRatePerPosTime=data.rightLocalGaussSpikeRatePerPosTime/GAUSSKERNELWIDTHSEC; %~spikes per sigma width of gaussian kernel (Hz)
        else
            %get left rate in different way, interp from original saved in
            %directionSpecificStats
            leftGaussTimeAxisTEMP=data.directionSpecificStats.(currFieldDirection).gaussTimeAxis;
            leftLocalGaussSpikeRatePerTimeTEMP=data.directionSpecificStats.(currFieldDirection).localGaussSpikeRatePerTime;
            localGaussSpikeRatePerPosTime=interp1(leftGaussTimeAxisTEMP,leftLocalGaussSpikeRatePerTimeTEMP,posInfo.positionTimeAxisSec);
            
            %localGaussSpikeRatePerPosTime=data.leftLocalGaussSpikeRatePerPosTime/GAUSSKERNELWIDTHSEC; %~spikes per sigma width of gaussian kernel (Hz)
        end
        
        if(max(localGaussSpikeRatePerPosTime)==0 || length(localGaussSpikeRatePerPosTime)==1)
               localGaussSpikeRatePerPosTime=NaN(size(localNormTimeInFieldPerTime));
        end
        %localGaussSpikeRatePerPosTime=data.localGaussSpikeRatePerPosTime;
        %else
        %    localGaussSpikeRatePerPosTime=NaN(size(localNormTimeInFieldPerTime));
        %end
        
        localSpikePhasePerTime=mod(localSpikePhasePerTime+offsetAngle,360);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %remove backwards moving spacetime points (another map)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        localSpikePhasePerTime(goingBackwardIdxes)=NaN;
        localGaussSpikeRatePerPosTime(goingBackwardIdxes)=NaN;
        
        localDistInFieldPerTime(goingBackwardIdxes)=NaN;
        localTimeInFieldPerTime(goingBackwardIdxes)=NaN;
        localNormDistInFieldPerTime(goingBackwardIdxes)=NaN;
        localNormTimeInFieldPerTime(goingBackwardIdxes)=NaN;
        
        absSpeedPerTimeStep=abs(speedPerTimeStep);

        %scatter(localNormDistInFieldPerTime,localNormTimeInFieldPerTime,5,localSpikePhasePerTime)
        %{
        if(length(localTimeInFieldPerTime) ~= length(localSpikePhasePerTime) || length(localTimeInFieldPerTime) ~= length(localGaussSpikeRatePerPosTime))
            close all
            continue
        end
        %}
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %circ shift phase s.t. starts at 360
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        if(contains(fileBaseName,'atsby'))
             bestShiftDeg=180;
             
        elseif(contains(fileBaseName,'chilles'))
            bestShiftDeg=-45; %45
            
        elseif(contains(fileBaseName,'uddy'))
            bestShiftDeg=-45;
            
         elseif(contains(fileBaseName,'icero'))
            bestShiftDeg=-45;
            
        else
             
             %bestShiftDeg=180;
              bestShiftDeg=0;
              %close all
              %continue
         end
        %first coarse then fine circ shift adjustment
        localSpikePhasePerTime=mod(localSpikePhasePerTime-bestShiftDeg,360);
                %figure(1); plot(localNormDistInFieldPerTime,localSpikePhasePerTime,'k.'); ylim([0 360]); xlim([0 1])

        
        max(localSpikePhasePerTime)
        
        localNormDistInFieldPerTime(:)=localNormDistInFieldPerTime(:)-zeroVal;
        localNormTimeInFieldPerTime(:)=localNormTimeInFieldPerTime(:)-zeroVal;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %TAKE OUT LOW SPEED SPIKES...
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         minSpeed=0.05;
    lowSpeedIdxes=abs(speedPerTimeStep(:))<minSpeed;
       localSpikePhasePerTime(lowSpeedIdxes)=NaN;
    localGaussSpikeRatePerPosTime(lowSpeedIdxes)=NaN;
    
        localSpikePhaseJumpPerTime=[NaN; angdiffDeg(localSpikePhasePerTime(:)) ];
        localSpikePhaseJumpPerTime(localSpikePhaseJumpPerTime==0)=NaN; %~same theta cycle time steps
      
        currNumFields=max(data.directionSpecificStats.(currFieldDirection).fieldNumPerTime);
        
        if(isnan(currNumFields))
            currNumFields=0;
            continue
        end
        fieldCount=fieldCount+currNumFields;
        cellCount=cellCount+1;
        
         fieldNumPerTime=data.directionSpecificStats.(currFieldDirection).fieldNumPerTime;
         
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %lap numbering comes from single reference unit!! 5/12/2021 TJ
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         lapNumPerTime=refUnitData.directionSpecificStats.(currFieldDirection).lapNumberPerTime;
         
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %get average speed per field lap and average phase per field lap
        %(log latency model predicts these related)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        avgFieldSpeedPerLap=abs(data.directionSpecificStats.(currFieldDirection).fieldWidthsPerLap./data.directionSpecificStats.(currFieldDirection).fieldDurationsPerLap);
      
       if(isempty(avgFieldSpeedPerLap))
           continue
       end
        %currNumLaps=length(avgFieldSpeedPerLap);
        currNumLaps=length(lapStartTimes);
        
        circAvgPhaseInFieldPerLap=NaN(currNumFields,currNumLaps);
        timeSinceExpStartPerLap=NaN(currNumFields,currNumLaps);
        
        avgTimeToFieldEndPerLap=NaN(currNumFields,currNumLaps);
        avgTimeInFieldPerLap=NaN(currNumFields,currNumLaps);
        circMedianPhaseInFieldPerLap=NaN(currNumFields,currNumLaps);
        phaseTimeLinFitSlopeInFieldPerLap=NaN(currNumFields,currNumLaps);
        phaseTimeLinFitOffsetInFieldPerLap=NaN(currNumFields,currNumLaps);
        isGoodLapForEachField=data.directionSpecificStats.(currFieldDirection).isGoodLapForEachField;
        
        masterPhaseJumpAmts=[];
        masterPhaseJumpInFieldTime=[];
        masterPhaseJumpInFieldDist=[];
        masterPhaseJumpSpeeds=[];
        masterPhaseJumpLapNum=[];
        masterPhaseJumpFieldDur=[];
        masterPhaseJumpStartPhase=[];
                jumpNoiseThresh=10; %deg
                minPrecessionJump=-45;
                 minPrecessionJump=-40;
                 minPrecessionJump=0;
                minSeqSpeed=0.35; %m/s
                 minSeqSpeed=0.1; %m/s
                 minSeqSpeed=0.05; %m/s
        %minSeqSpeed=0.2; %m/s
        
        %useRealTime=1;
        useRealTime=0;
        %{
        figure;
        showFieldImage(data,di)
        %}
        
        avgFieldSpeedPerLap=NaN(currNumFields,currNumLaps);
        for fii=1:currNumFields  
            for li=1:currNumLaps
                %{
                if(~isGoodLapForEachField(fii,li))
                    continue
                end
                %}
                %allTimeIdxesInCurrFieldAndLap=(fieldNumPerTime==fii) & (lapNumPerTime==li);
                timeIdxesInLap=posInfo.positionTimeAxisSec>=lapStartTimes(li) & posInfo.positionTimeAxisSec<=lapStopTimes(li);
                
                allTimeIdxesInCurrFieldAndLap=(fieldNumPerTime==fii) & (timeIdxesInLap);
                %{
                figure(190)
                plot(data.positionTimeAxis,timeIdxesInLap,'k')
                hold on
                plot(data.positionTimeAxis,allTimeIdxesInCurrFieldAndLap,'r')
                
                if(li==1)
                for sii=1:length(unitDirSpikeTimes)
                    plot([unitDirSpikeTimes(sii) unitDirSpikeTimes(sii)],ylim,'b')
                    hold on
                end
                end
                %}
                
                currFieldLapTimes=posInfo.positionTimeAxisSec(allTimeIdxesInCurrFieldAndLap);
                currSessionExpStartTime=posInfo.positionTimeAxisSec(1);
                
                if(isempty(currFieldLapTimes))
                    continue
                end
                currFieldLapEntryTime=min(currFieldLapTimes);
                currFieldLapExitTime=max(currFieldLapTimes);
                
                %find speed range (variability) this lap to restrict to
                %straight speed laps
               %speedRangeThisFieldLap=data.thetaBasedTimeVsPhaseInfo.(currFieldDirection).(sprintf('field%d',fii)).allLapInFieldSpeedRangePerTraversal(li);
                currFieldLapSpeeds=abs(posInfo.speedAlongTrackPerTimeMSec(allTimeIdxesInCurrFieldAndLap));
                speedRangeThisFieldLap=range(currFieldLapSpeeds);
                            
                  
                expectedFieldStart=data.(sprintf('%sFieldStartEndM',currFieldDirection))(fii);
                expectedFieldEnd=data.(sprintf('%sFieldStartEndM',currFieldDirection))(fii+1);
                
                expectedFieldWidth=abs(expectedFieldEnd-expectedFieldStart);
                
                %avgFieldSpeedPerLap(fii,li)=expectedFieldWidth/(currFieldLapExitTime-currFieldLapEntryTime);
                
                avgFieldSpeedPerLap(fii,li)=nanmean(abs(speedPerTimeStep(allTimeIdxesInCurrFieldAndLap)));
                
                
                
                
                if(speedRangeThisFieldLap>maxSpeedRangeInField)
                    localSpikePhasePerTime(allTimeIdxesInCurrFieldAndLap)=NaN;
                    localGaussSpikeRatePerPosTime(allTimeIdxesInCurrFieldAndLap)=NaN;
                    continue %to next lap
                end
                
                hasMinSeqSpeedIdxesPerTime=abs(speedPerTimeStep(:)')>minSeqSpeed;
                
                 if(length(fieldNumPerTime) ~= length(speedPerTimeStep(:)))
                    continue
                 end
                 
                
                
                 currFieldLapDirSpikeTimeIdxes=unitDirSpikeTimes>=currFieldLapEntryTime & unitDirSpikeTimes<=currFieldLapExitTime;
                 
                 currFieldLapDirSpikeTimes=unitDirSpikeTimes(unitDirSpikeTimes>=currFieldLapEntryTime & unitDirSpikeTimes<=currFieldLapExitTime);
                currFieldLapDirSpikeRefPhases=unitDirSpikeRefPhases(currFieldLapDirSpikeTimeIdxes);
                currFieldLapDirSpikeCycleIDs=unitDirSpikeRefCycleIDs(currFieldLapDirSpikeTimeIdxes);
                
                 currFieldLapSpikeTimesFromFieldStart=currFieldLapDirSpikeTimes-currFieldLapEntryTime;
                 currFieldLapSpikeTimesToFieldEnd=currFieldLapExitTime-currFieldLapDirSpikeTimes;
                 
                 currFieldLapSpikeTimesFromExpStart=currFieldLapDirSpikeTimes-currSessionExpStartTime;
                
                %{
                figure(191)
                plot(unitDirSpikeTimes,'ko'); hold on; plot(xlim,[currFieldLapEntryTime currFieldLapEntryTime]);  plot(xlim,[currFieldLapExitTime currFieldLapExitTime])
                %}
                uniqueCycleIDs=unique(currFieldLapDirSpikeCycleIDs);
                uniqueCycleAvgPhases=NaN(length(uniqueCycleIDs),1);
                uniqueCycleAvgTimesFromStart=NaN(length(uniqueCycleIDs),1);
                uniqueCycleAvgTimesFromExpStart=NaN(length(uniqueCycleIDs),1);
                uniqueCycleAvgTimesToEnd=NaN(length(uniqueCycleIDs),1);
                for uci=1:length(uniqueCycleIDs)
                    currCycleID=uniqueCycleIDs(uci);
                    
                    currCycleIDspikeIdxes=currFieldLapDirSpikeCycleIDs==currCycleID;
                    uniqueCycleAvgPhases(uci)=nanmean(currFieldLapDirSpikeRefPhases(currCycleIDspikeIdxes));
                    uniqueCycleAvgTimesFromStart(uci)=nanmean(currFieldLapSpikeTimesFromFieldStart(currCycleIDspikeIdxes));
                    uniqueCycleAvgTimesToEnd(uci)=nanmean(currFieldLapSpikeTimesToFieldEnd(currCycleIDspikeIdxes));
                    uniqueCycleAvgTimesFromExpStart(uci)=nanmean(currFieldLapSpikeTimesFromExpStart(currCycleIDspikeIdxes));
                end
                
                %{
                if(length(currFieldLapDirSpikeCycleIDs)>5)
                    disp('')
                    
                else
                    continue
                end
                %}
                
                
                
                allSpikePhasesInCurrFieldAndLap=uniqueCycleAvgPhases;
                allSpikeExpTimesInCurrFieldAndLap=uniqueCycleAvgTimesFromExpStart;
                
                %allSpikePhasesInCurrFieldAndLap=localSpikePhasePerTime(allTimeIdxesInCurrFieldAndLap);
                 allSpikePhaseJumpsInCurrFieldAndLap=localSpikePhaseJumpPerTime(allTimeIdxesInCurrFieldAndLap & hasMinSeqSpeedIdxesPerTime);
                
                %timesInCurrFieldAndLap= data.positionTimeAxis(allTimeIdxesInCurrFieldAndLap);
                if(useRealTime)
                    timesInCurrFieldAndLap= localTimeInFieldPerTime(allTimeIdxesInCurrFieldAndLap);
                else
                    timesInCurrFieldAndLap= localNormTimeInFieldPerTime(allTimeIdxesInCurrFieldAndLap);
                end

                distsInCurrFieldAndLap=localNormDistInFieldPerTime(allTimeIdxesInCurrFieldAndLap);
                
                speedsInCurrFieldAndLap=absSpeedPerTimeStep(allTimeIdxesInCurrFieldAndLap);
                startPhasesInCurrFieldAndLap=localSpikePhasePerTime(allTimeIdxesInCurrFieldAndLap);
                
                notNaNIdxesPhase=~isnan(allSpikePhasesInCurrFieldAndLap);
                
                if(~isempty(allSpikePhasesInCurrFieldAndLap(notNaNIdxesPhase)))
                    %{
                    linCoeffs=polyfit(timesInCurrFieldAndLap(notNaNIdxesPhase),allSpikePhasesInCurrFieldAndLap(notNaNIdxesPhase),1);
                    inFieldLapSlope=linCoeffs(1);
                    inFieldLapOffset=linCoeffs(2);
                    phaseTimeLinFitSlopeInFieldPerLap(fii,li)=inFieldLapSlope;
                    phaseTimeLinFitOffsetInFieldPerLap(fii,li)=inFieldLapOffset;
                    %}
                    
                    circAvgPhaseInFieldPerLap(fii,li)=circMeanDeg(allSpikePhasesInCurrFieldAndLap(:));
                    timeSinceExpStartPerLap(fii,li)=nanmean(allSpikeExpTimesInCurrFieldAndLap);
                    
                    avgTimeInFieldPerLap(fii,li)=nanmean(uniqueCycleAvgTimesFromStart(:));
                    avgTimeToFieldEndPerLap(fii,li)=nanmean(uniqueCycleAvgTimesToEnd(:));
                    
                    circMedianPhaseInFieldPerLap(fii,li)=circMedianDeg(allSpikePhasesInCurrFieldAndLap(:));
                 %continue
                end
               
                %notNaNIdxesJumps=~isnan(allSpikePhaseJumpsInCurrFieldAndLap) & allSpikePhaseJumpsInCurrFieldAndLap<-jumpNoiseThresh & allSpikePhaseJumpsInCurrFieldAndLap>(-180+jumpNoiseThresh); 
                %notNaNIdxesJumps=~isnan(allSpikePhaseJumpsInCurrFieldAndLap) & allSpikePhaseJumpsInCurrFieldAndLap<-(jumpNoiseThresh+30) & allSpikePhaseJumpsInCurrFieldAndLap>(-180+jumpNoiseThresh); 
                 notNaNIdxesJumps=~isnan(allSpikePhaseJumpsInCurrFieldAndLap) & allSpikePhaseJumpsInCurrFieldAndLap<minPrecessionJump; 
     
                 currFieldDur=(data.directionSpecificStats.(currFieldDirection).avgFieldDurations(fii));
                %drastic and meager shifts seem contaminated with noise
                %min speed 40cm/s
                
                nonNaNPhaseJumps=allSpikePhaseJumpsInCurrFieldAndLap(notNaNIdxesJumps);
                nonNaNPhaseJumpTimes=timesInCurrFieldAndLap(notNaNIdxesJumps);
                nonNaNPhaseJumpDists=distsInCurrFieldAndLap(notNaNIdxesJumps);
                nonNaNPhaseJumpSpeeds=speedsInCurrFieldAndLap(notNaNIdxesJumps);
                nonNaNPhaseJumpStartPhases=startPhasesInCurrFieldAndLap(notNaNIdxesJumps);
                
                if(~isempty(nonNaNPhaseJumps))
                    masterPhaseJumpAmts=[masterPhaseJumpAmts; nonNaNPhaseJumps(:)];
                    masterPhaseJumpInFieldTime=[masterPhaseJumpInFieldTime; nonNaNPhaseJumpTimes(:)];
                    masterPhaseJumpInFieldDist=[masterPhaseJumpInFieldDist; nonNaNPhaseJumpDists(:)];
                    masterPhaseJumpLapNum=[masterPhaseJumpLapNum; repelem(li,length(nonNaNPhaseJumps(:)))'];
                    masterPhaseJumpFieldDur=[masterPhaseJumpFieldDur; repelem(currFieldDur,length(nonNaNPhaseJumps(:)))'];
                    masterPhaseJumpSpeeds=[masterPhaseJumpSpeeds; nonNaNPhaseJumpSpeeds(:)];
                    masterPhaseJumpStartPhase=[masterPhaseJumpStartPhase; nonNaNPhaseJumpStartPhases(:)];
                end    
            end
        end
        
        phaseJumpVsTimeVsDistInField=[masterPhaseJumpInFieldDist(:),masterPhaseJumpInFieldTime(:), masterPhaseJumpAmts(:),masterPhaseJumpLapNum(:),masterPhaseJumpFieldDur(:),masterPhaseJumpSpeeds(:),masterPhaseJumpStartPhase(:)];

        if(size(phaseJumpVsTimeVsDistInField,1)<2) 
            continue
        end
        
        if(useRealTime)         
           corrIdxes=~isnan(phaseJumpVsTimeVsDistInField(:,2)) & ~isnan(phaseJumpVsTimeVsDistInField(:,3)) & phaseJumpVsTimeVsDistInField(:,2)<=phaseJumpVsTimeVsDistInField(:,5); %only within typical field duration
        else
         corrIdxes=~isnan(phaseJumpVsTimeVsDistInField(:,2)) & ~isnan(phaseJumpVsTimeVsDistInField(:,3)) & phaseJumpVsTimeVsDistInField(:,2)<=1; %only within typical field duration
        end
        
        if(length(phaseJumpVsTimeVsDistInField(corrIdxes,2))<2)
            continue
        end
        rMat=corrcoef([phaseJumpVsTimeVsDistInField(corrIdxes,2),phaseJumpVsTimeVsDistInField(corrIdxes,3)]);
        
        rJumpTime=rMat(1,2);
        
        linCoeffsJump=polyfit(phaseJumpVsTimeVsDistInField(corrIdxes,2),phaseJumpVsTimeVsDistInField(corrIdxes,3),1);
                    jumpVsTimeSlope=linCoeffsJump(1);
                    jumpVsTimeOffset=linCoeffsJump(2);
                    
        if(showPlots)
        figure(1010);
        %{
        subplot(1,2,2)
        plot(phaseJumpVsTimeVsDistInField(:,1),phaseJumpVsTimeVsDistInField(:,3),'k.','MarkerSize',10)
        xlim([0 1.1])
               ylim([-180 0])
        title('phase jump  cycle vs distance in field')
             subplot(1,2,1)
         %}
        
       subplot(2,1,1)
        
        plot(phaseJumpVsTimeVsDistInField(corrIdxes,2),phaseJumpVsTimeVsDistInField(corrIdxes,3),'k.','MarkerSize',30)
                %scatter(phaseJumpVsTimeVsDistInField(:,2),phaseJumpVsTimeVsDistInField(:,3),500,phaseJumpVsTimeVsDistInField(:,4),'filled')
                %colormap(gca,jet)
                %cbJ=colorbar;
                %ylabel(cbJ,'lap no.')

                
                if(useRealTime)
                    xlabel('Time in field (sec)')
                else
                    xlabel('Time in field (frac)')
                end
                ylabel('Phase jump from last cycle (deg)')
            xlim([0 1])
            %ylim([-180 0])
            ylim([-180 minPrecessionJump])
            
            %daspect([1 360 1])
            if(useRealTime)
                slopeUnitStr='deg-change/cycle/sec';
            else
                slopeUnitStr='deg-change/cycle/field';
            end
        title({removeUnderscores(fileBaseName),'cycle-to-cycle phase jump vs time in field',sprintf('R=%.3f, slope=%.2f %s , min speed: %.2f m/s',rJumpTime,jumpVsTimeSlope,slopeUnitStr,minSeqSpeed)})
       
        maxFigHalfWidth
   
        
         %uberTitle(removeUnderscores(fileBaseName))
              setFigFontTo(16)
              
              
        end
        
        speedVsMeanPhaseStats=[]; %reset variable between loops, continues
        speedVsMeanPhaseStats.(currFieldDirection).circAvgPhaseInFieldPerLap=circAvgPhaseInFieldPerLap;
        speedVsMeanPhaseStats.(currFieldDirection).avgTimeInFieldPerLap=avgTimeInFieldPerLap;
        speedVsMeanPhaseStats.(currFieldDirection).avgTimeToFieldEndPerLap=avgTimeToFieldEndPerLap;
        
        speedVsMeanPhaseStats.(currFieldDirection).circMedianPhaseInFieldPerLap=circMedianPhaseInFieldPerLap;
        speedVsMeanPhaseStats.(currFieldDirection).timeSinceExpStartPerLap=timeSinceExpStartPerLap;
        %{
         speedVsMeanPhaseStats.(currFieldDirection).phaseTimeLinFitSlopeInFieldPerLap=phaseTimeLinFitSlopeInFieldPerLap;
          speedVsMeanPhaseStats.(currFieldDirection).phaseTimeLinFitOffsetInFieldPerLap=phaseTimeLinFitOffsetInFieldPerLap;
          %}
                speedVsMeanPhaseStats.(currFieldDirection).currNumFields=currNumFields;
        speedVsMeanPhaseStats.(currFieldDirection).currNumLaps=currNumLaps;
        speedVsMeanPhaseStats.(currFieldDirection).avgFieldSpeedPerLap=avgFieldSpeedPerLap;
        speedVsMeanPhaseStats.(currFieldDirection).fieldWidthsPerLap=data.directionSpecificStats.(currFieldDirection).fieldWidthsPerLap;
        speedVsMeanPhaseStats.(currFieldDirection).fieldDurationsPerLap=data.directionSpecificStats.(currFieldDirection).fieldDurationsPerLap;
        
        %if(autoShiftOffset==0)
            
            if(di==1 && strcmp(currFieldDirection,'rightward'))
                speedVsMeanPhaseStatsRightwardNoAutoShiftRef=speedVsMeanPhaseStats;
                %save(currFilePath,'speedVsMeanPhaseStatsRightwardNoAutoShift','-append')
                save(currFilePath,'speedVsMeanPhaseStatsRightwardNoAutoShiftRef','-append')
            elseif(di==2 && strcmp(currFieldDirection,'leftward'))
                speedVsMeanPhaseStatsLeftwardNoAutoShiftRef=speedVsMeanPhaseStats;
                  %save(currFilePath,'speedVsMeanPhaseStatsLeftwardNoAutoShift','-append')
                  save(currFilePath,'speedVsMeanPhaseStatsLeftwardNoAutoShiftRef','-append')
            end
            %close all
            %continue
        %else
             if(di==1 && strcmp(currFieldDirection,'rightward'))
                phaseJumpVsTimeVsDistInFieldRightward=phaseJumpVsTimeVsDistInField;
                rPhaseJumpVsTimeRight=rJumpTime;
                jumpVsTimeSlopeRight=jumpVsTimeSlope;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %RESET
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %{
                phaseJumpVsTimeVsDistInFieldRightward=NaN;
                rPhaseJumpVsTimeRight=NaN;
                jumpVsTimeSlopeRight=NaN;
                %}
                
                save(currFilePath,'phaseJumpVsTimeVsDistInFieldRightward','rPhaseJumpVsTimeRight','jumpVsTimeSlopeRight','minSeqSpeed','-append')
              
                %save(currFilePath,'phaseJumpVsTimeVsDistInFieldRightward','rPhaseJumpVsTimeRight','jumpVsTimeSlopeRight','minSeqSpeed','-append')

            elseif(di==2 && strcmp(currFieldDirection,'leftward'))
                phaseJumpVsTimeVsDistInFieldLeftward=phaseJumpVsTimeVsDistInField;
                rPhaseJumpVsTimeLeft=rJumpTime;
                jumpVsTimeSlopeLeft=jumpVsTimeSlope;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %RESET
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %{
                phaseJumpVsTimeVsDistInFieldLeftward=NaN;
                rPhaseJumpVsTimeLeft=NaN;
                jumpVsTimeSlopeLeft=NaN;
                %}
                  save(currFilePath,'phaseJumpVsTimeVsDistInFieldLeftward','rPhaseJumpVsTimeLeft','jumpVsTimeSlopeLeft','minSeqSpeed','-append')
            end
        %end
        %continue
        
        
       if(autoShiftOffset)
         [bestShiftDeg]=getBestShiftDeg(localSpikePhasePerTime,localNormDistInFieldPerTime,localNormTimeInFieldPerTime)
        localSpikePhasePerTime=mod(localSpikePhasePerTime-bestShiftDeg,360);
       end
        
        fJ=figure(1010);
        %close(fJ)
        %continue
        
        %fT=figure(1011); 
        phaseChanges=[NaN diff(localSpikePhasePerTime)];
        changeIdxes=phaseChanges~=0;
        subplot(2,1,2)
        plot(localNormTimeInFieldPerTime(changeIdxes),localSpikePhasePerTime(changeIdxes),'k.','MarkerSize',25);ylim([0 360]); xlim([0 1])
        xlabel('Time in field (frac)')
        ylabel('Spike phase (degrees)')
        maxFigHalfWidth
        setFigFontTo(18)
        title(removeUnderscores(fileBaseName))
        %daspect([1 360*1.5 1])
        
        saveas(gcf,fullfile(savePhaseJumpVsTimeDir,sprintf('%s_PhaseJumpVsTime_%s.tif',fileBaseName,currFieldDirection)))

        %saveas(gcf,fullfile(savePhaseJumpVsTimeDir,sprintf('%s_examplePhaseVsTimeSingleCell_%s.tif',fileBaseName,currFieldDirection)))
        
        %close(fT)
        
        masterSpikePhasePerTime=[masterSpikePhasePerTime;localSpikePhasePerTime(:)];
         masterDistPerTime=[masterDistPerTime; localNormDistInFieldPerTime(:)];
    masterLapTimePerTime=[masterLapTimePerTime; localNormTimeInFieldPerTime(:)];
    masterSpikeRatePerTime=[masterSpikeRatePerTime;localGaussSpikeRatePerPosTime(:)];
    %masterSpeedPerTime=[masterSpeedPerTime; abs(data.speedPerTimeStep(:))];
        masterSpeedPerTime=[masterSpeedPerTime; absSpeedPerTimeStep(:)];
        
        totalUsedFieldCount=totalUsedFieldCount+1;
        
         if(length(fieldNumPerTime)~=length(localSpikePhasePerTime(:)))
             fieldNumPerTime=NaN(size(localSpikePhasePerTime));
         end
        masterFieldNumPerTime=[masterFieldNumPerTime; fieldNumPerTime(:)+fieldCount-1];
        
       % masterSesNumPerTime=[masterSesNumPerTime; repmat(currSessionNum,size(localSpikePhasePerTime(:)))];
        
        masterDirectionPerTime=[masterDirectionPerTime; repmat(di,size(localSpikePhasePerTime(:)))];
    
         %plotSpacetimePhaseHeatMap([localNormDistInFieldPerTime(:), localNormTimeInFieldPerTime(:), localSpikePhasePerTime(:)],0.01,0);
         %plotSpacetimePhaseHeatMap([localNormDistInFieldPerTime(:), localNormTimeInFieldPerTime(:), localSpikePhasePerTime(:)],0.02,0);
        max(localGaussSpikeRatePerPosTime)
         plotSpacetimePhaseHeatMap([localNormDistInFieldPerTime(:), localNormTimeInFieldPerTime(:), localSpikePhasePerTime(:)],heatResUnit,0,[],localGaussSpikeRatePerPosTime(:),totalUsedFieldCount);

         %scatter(localDistInFieldPerTime,localTimeInFieldPerTime,20,localSpikePhasePerTime,'filled')
        %scatter(localNormDistInFieldPerTime,localNormTimeInFieldPerTime,20,360-localSpikePhasePerTime,'filled')
        hold on
        colormap(jet)
        colorbar
        
        %xlabel('Distance in field (m)')
        %ylabel('Time in field (sec)')
        xlabel('Distance in field (frac of avg width)')
        ylabel('Time in field (frac of avg duration)')
        uberTitle(sprintf('%s, %s',removeUnderscores(fileBaseName),currFieldDirection))
        setFigFontTo(18)
        
      
        
        %{
        avgFieldWidths=data.directionSpecificStats.(currFieldDirection).avgFieldWidths;
        avgFieldDurations=data.directionSpecificStats.(currFieldDirection).avgFieldDurations;
        for i=1:length(avgFieldWidths)
            plot([avgFieldWidths(i) avgFieldWidths(i)],ylim,'k--','LineWidth',3)
            plot(xlim,[avgFieldDurations(i) avgFieldDurations(i)],'k--','LineWidth',3)
        end
        
        %}
        %xlim([0 1.1])
        %ylim([0 1.5])
        %{
            xlim([0 1.1])
            ylim([0 1.1])
          plot(xlim,[1 1],'k--','LineWidth',3)
         plot([1 1],ylim, 'k--','LineWidth',3)
        %}
        subplot(1,2,1)
        caxis([0 360])
        daspect([1 1 1])
        %plotHorzLinesEvery05
        set(gca,'XMinorTick','on','YMinorTick','on')
       
        
        saveas(gcf,fullfile(spaceTimeResponsesDir,sprintf('%s_SpaceTimeResponses_%s.tif',fileBaseName,currFieldDirection)))
         
        %{
        if(iscell( data.imagePaths))
            for ii=1:length(data.imagePaths)
                figure(ii+1)
                imshow(data.imagePaths{ii})
            end
        else
            figure(2)
                imshow(data.imagePaths)
        end
        %}
       autoArrangeFigures()
        close all
    end
end
masterDistPerTime(:)=masterDistPerTime(:);
masterLapTimePerTime(:)=masterLapTimePerTime(:);
%[latestOutputDataStruct]=plotSpacetimePhaseHeatMap([masterDistPerTime(:), masterLapTimePerTime(:), masterSpikePhasePerTime(:)],0.01,1);
[latestOutputDataStruct]=plotSpacetimePhaseHeatMap([masterDistPerTime(:), masterLapTimePerTime(:), masterSpikePhasePerTime(:)],heatRes,1,[],masterSpikeRatePerTime(:),totalUsedFieldCount);

save('masterSpacetimePhaseResponse','-v7.3')
else
    load('masterSpacetimePhaseResponse.mat')
    
    if(~exist('latestOutputDataStruct','var'))
        [latestOutputDataStruct]=plotSpacetimePhaseHeatMap([masterDistPerTime(:), masterLapTimePerTime(:), masterSpikePhasePerTime(:)],heatRes,1,[],masterSpikeRatePerTime(:),totalUsedFieldCount);
    end
    
    
    
    bothZeroIdxes=masterDistPerTime==0 & masterLapTimePerTime==0;
    noSpikesIdxes=isnan(masterSpikePhasePerTime);
    
    masterDistPerTime(bothZeroIdxes)=NaN;
    masterLapTimePerTime(bothZeroIdxes)=NaN;
    
     masterDistPerTime(noSpikesIdxes)=NaN;
    masterLapTimePerTime(noSpikesIdxes)=NaN;
    
   
    
    distTimesTime=masterDistPerTime.*masterLapTimePerTime;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %plot joint probability slices holding distance or time into field 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x=masterDistPerTime;
    y=masterLapTimePerTime;
    z=masterSpikePhasePerTime/360;
    
    %plot3DjointDistrSlices(x,y,z,'Distance','Time','Phase')
     % plot3DjointDistrSlices(y,x,z,'Time','Distance','Phase')
     
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %analyze symmetric log shape along x=t line i.e. x-t=0 line
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    distBinsPhase=latestOutputDataStruct.distBinsPhase;
    timeBinsPhase=latestOutputDataStruct.timeBinsPhase;
    meanPhaseSmooth=latestOutputDataStruct.meanPhaseSmooth;
    dataCountPerBinPhase=latestOutputDataStruct.dataCountPerBinPhase;
    
    %meanPhaseSmooth=nanGaussSmooth(meanPhaseSmooth);
        meanPhaseSmooth=nanGaussSmoothCentral(meanPhaseSmooth);
    meanPhaseSmooth(dataCountPerBinPhase<10)=NaN;
    
       
    [X,T]=meshgrid(distBinsPhase,timeBinsPhase);
    centralSpeedLineOffsetUpper=0.2;
    centralSpeedLineOffsetLower=0.2;
    
    centralSpeedIdxes=X-T<centralSpeedLineOffsetUpper & X-T>-centralSpeedLineOffsetLower;
   
    centralSpeedPhases=nanGaussSmoothCentral(meanPhaseSmooth);
    centralSpeedPhases(~centralSpeedIdxes)=NaN;
    
    
    centralSpeedPhases=imrotate(centralSpeedPhases,45);
    centralSpeedPhases(centralSpeedPhases==0)=NaN;
    edgeCutoff=10;
        edgeCutoff=0;
    centralSpeedPhases(:,1:edgeCutoff)=NaN;
    centralSpeedPhases(:,(end-edgeCutoff+1):end)=NaN;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %plot central speed phases
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fC=figure
    subplot(3,1,1)
    %omarPcolor(distBinsPhase,timeBinsPhase,meanPhaseSmooth',fC)
        surf(distBinsPhase,timeBinsPhase,meanPhaseSmooth'/360)
    colormap(jet)
    colorbar
    daspect([1 1 1])
 
        
    hold on
    plot([0 1],[0 1],'k--','LineWidth',3)
    
      subplot(3,1,2)
    omarPcolor(1:length(centralSpeedPhases),1:length(centralSpeedPhases),centralSpeedPhases/360',fC)
    %imagesc(centralSpeedPhases/360)
    colormap(jet)
    colorbar
    daspect([1 1 1])
    
    hold on
    plot([0 1],[0 1],'k--','LineWidth',3)
    
    centralSpeedAvgPhase=nanmedian(centralSpeedPhases,1);
    

    subplot(3,1,3); plot(smooth(centralSpeedAvgPhase,7),'k-')
    xlim([0 Inf])
    axis square
    
    disp('')
    
    %distTimesTime=(masterDistPerTime+masterLapTimePerTime)/2;
    %{
    fH=figure; 
    subplot(1,3,1); %plot(distTimesTime,masterSpikePhasePerTime/360,'k.')
    
    binWidth=0.02;
    distTimesTime(distTimesTime>1)=NaN;
    plotJointHeatMap(distTimesTime,masterSpikePhasePerTime/360,binWidth,fH);
 
    xlim([0 1])
      ylim([0 1])
      caxis([0 150])
            daspect([1 1 1])
            plot([0 1],[1 0],'k--','LineWidth',3)
     subplot(1,3,2); %plot(masterDistPerTime,masterSpikePhasePerTime/360,'k.')
 
    masterDistPerTimeDisp=masterDistPerTime;
    masterDistPerTimeDisp(masterDistPerTimeDisp>1)=NaN;
    plotJointHeatMap(masterDistPerTimeDisp.^2,masterSpikePhasePerTime/360,binWidth,fH);
    xlim([0 1])
      ylim([0 1])
         caxis([0 150])
            daspect([1 1 1])
            plot([0 1],[1 0],'k--','LineWidth',3)
            

     subplot(1,3,3); %plot(masterLapTimePerTime,masterSpikePhasePerTime/360,'k.')

    masterLapTimePerTimeDisp=masterLapTimePerTime;
    masterDistPerTimeDisp(masterLapTimePerTime>1)=NaN;
    plotJointHeatMap(masterDistPerTimeDisp.^2,masterSpikePhasePerTime/360,binWidth,fH);
    xlim([0 1])
      ylim([0 1])
         caxis([0 150])
      
      daspect([1 1 1])

         plot([0 1],[1 0],'k--','LineWidth',3)
%}     

         %speedCutOff=1; %m/s %speedCutOff=prctile(x,99.5); %spuriously fast... 
     %x=x/speedCutOff;
     
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %time vs distance, holding phase
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     z=masterSpikePhasePerTime(:)/360;
     zStr=sprintf('Phase (frac of cycle)');
    y=masterDistPerTime(:);
      yStr=sprintf('Distance (frac of field)');
    x=masterLapTimePerTime(:);
        xStr=sprintf('Time (frac of field)');
        
        %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %time vs distance, holding speed
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     z=masterSpeedPerTime(:);
     zStr=sprintf('Speed (m/s)');
    y=masterDistPerTime(:);
    yStr=sprintf('Distance (frac of field)');
    x=masterLapTimePerTime(:);
        xStr=sprintf('Time (frac of field)');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %speed vs distance, holding time
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     x=masterSpeedPerTime(:);
     xStr=sprintf('Speed (m/s)');
    y=masterDistPerTime(:);
    yStr=sprintf('Distance (frac of field)');
    z=masterLapTimePerTime(:);
        zStr=sprintf('Time (frac of field)');
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %speed vs distance, holding phase
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      z=masterSpikePhasePerTime(:)/360;
     zStr=sprintf('Phase (frac of cycle)');
    y=masterDistPerTime(:);
    yStr=sprintf('Distance (frac of field)');
     x=masterSpeedPerTime(:); 
     xStr=sprintf('Speed (m/s)');
     %}
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %speed vs phase, holding distance
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      y=masterSpikePhasePerTime(:)/360;
     yStr=sprintf('Phase (frac of cycle)');
    z=masterDistPerTime(:);
    zStr=sprintf('Distance (frac of field)');
     x=masterSpeedPerTime(:); 
     xStr=sprintf('Speed (m/s)');
    
     %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %speed vs phase, holding time
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      y=masterSpikePhasePerTime(:)/360;
     yStr=sprintf('Phase (frac of cycle)');
    z=masterLapTimePerTime(:);
        zStr=sprintf('Time (frac of field)');
     x=masterSpeedPerTime(:); 
     xStr=sprintf('Speed (m/s)');
     %}
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %speed vs phase, holding distance
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      y=masterSpikePhasePerTime(:)/360;
     yStr=sprintf('Phase (frac of cycle)');
    z=masterDistPerTime(:);
    zStr=sprintf('Distance (frac of field)');
    %z(z>0.3333)=NaN;
     x=masterSpeedPerTime(:); 
     xStr=sprintf('Speed (m/s)');
     
     %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %speed vs time, holding distance
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     x=masterSpeedPerTime(:);
     xStr=sprintf('Speed (m/s)');
    z=masterDistPerTime(:);
    zStr=sprintf('Distance (frac of field)');
    y=masterLapTimePerTime(:);
        yStr=sprintf('Time (frac of field)');
     %}
     
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %time vs phase, holding distance
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    z=masterDistPerTime(:);
       zStr=sprintf('Distance (frac of field)');
         x=masterLapTimePerTime(:);
        xStr=sprintf('Time (frac of field)');
            y=masterSpikePhasePerTime(:)/360;
    yStr=sprintf('Phase (frac of cycle)');
    
    
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %distance vs phase, holding time
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x=masterDistPerTime(:);
       xStr=sprintf('Distance (frac of field)');
         z=masterLapTimePerTime(:);
        zStr=sprintf('Time (frac of field)');
            y=masterSpikePhasePerTime(:)/360;
    yStr=sprintf('Phase (frac of cycle)');
    
    timeOffset=0.05;
    timeOffset=0.2;
    %masterLogTimeToEndPerTime=getLogTime(masterLapTimePerTime(:) - timeOffset);
    
    masterLogTimeToEndPerTime=masterLapTimePerTime(:);
    
  % getConditionalMutualInfo(x,y,z,xStr,yStr,zStr)
    
   % getLinearXYCoeffsforZstrata(x,y,z)
    %plotSpacetimePhaseHeatMap([masterDistPerTime(:), masterLapTimePerTime(:), masterSpikePhasePerTime(:)],0.01,1);
    %plotSpacetimePhaseHeatMap([masterDistPerTime(:), masterLapTimePerTime(:), masterSpikePhasePerTime(:)],0.01,1,latestOutputDataStruct);
        %plotSpacetimePhaseHeatMap([masterDistPerTime(:), masterLapTimePerTime(:), masterSpikePhasePerTime(:)],heatRes,1,latestOutputDataStruct,masterSpikeRatePerTime(:));
        %plotSpacetimePhaseHeatMap([masterDistPerTime(:), masterLapTimePerTime(:), masterSpikePhasePerTime(:)],heatRes,1,latestOutputDataStruct,masterSpikeRatePerTime(:));
        %heatRes=0.1;
        heatRes=0.01;
        plotSpacetimePhaseHeatMap([masterDistPerTime(:), masterLogTimeToEndPerTime(:), masterSpikePhasePerTime(:)],heatRes,1,[],masterSpikeRatePerTime(:),totalUsedFieldCount);

        
        
        %plotSpacetimePhaseHeatMap([masterDistPerTime(:), masterLapTimePerTime(:), masterSpikePhasePerTime(:)],0.025,1,latestOutputDataStruct);

    %scatter(localDistInFieldPerTime,localTimeInFieldPerTime,20,localSpikePhasePerTime,'filled')
        %scatter(localNormDistInFieldPerTime,localNormTimeInFieldPerTime,20,360-localSpikePhasePerTime,'filled')
        hold on
        colormap(jet)
        colorbar
        
        %xlabel('Distance in field (m)')
        %ylabel('Time in field (sec)')
        xlabel('Distance in field (frac of avg width)')
        ylabel('Time in field (frac of avg duration)')
        
        setFigFontTo(18)
        
        plot(xlim,[1 1],'k--','LineWidth',3)
         plot([1 1],ylim, 'k--','LineWidth',3)
         uberTitle('Spacetime field average summary')
         xlim([0 1.1])
         ylim([0 1.1])
         daspect([1 1 1])
         maxFig
         saveas(gcf,sprintf('HC3_SummaryAvgSpaceTimeResponse.tif'))
         %%
         timeFig=figure(111); 
         distFig=figure(1111); 
         logFitFigTime=figure(1112);
         logFitFigDist=figure(1113);
         %plot(masterDistPerTime,masterSpikePhasePerTime,'k.')
         numSpeedStrata=4;
         %numSpeedStrata=9;
          %numSpeedStrata=16;
           %   numSpeedStrata=8;
            numSpeedStrata=1;
         
         speedEdges=linspace(0,1.2,numSpeedStrata+1);
         
         speedEdges=linspace(0.05,1.2,numSpeedStrata+1);
          %speedEdges=linspace(0.5,1.2,numSpeedStrata+1);
         %speedEdges=linspace(0.2,1.2,numSpeedStrata+1);
         %speedEdges=linspace(0.35,1.2,numSpeedStrata+1);
         numFieldsAnalyzed=max(masterFieldNumPerTime);
   
             
         for speedStratum=1:numSpeedStrata
             rhoCircLinAllFields=[];
             rhoCircLogAllFields=[];
             allAsymIdx=[];
               %for currFieldNum=1:numFieldsAnalyzed
         masterDistPerTime(masterDistPerTime>(1-zeroVal))=NaN;
         
         %maxSpeed=1.2;
          %   minSpeed=0.05;
          maxSpeed=speedEdges(speedStratum+1);
             minSpeed=speedEdges(speedStratum);
             
             %{
         highSpeedThresh=0.9;
         highSpeedIdxes=masterSpeedPerTime<maxSpeed & masterSpeedPerTime>highSpeedThresh;
         lowSpeedThresh=0.5;
          lowSpeedIdxes=masterSpeedPerTime<lowSpeedThresh & masterSpeedPerTime>minSpeed;
         onlyHighSpeed=0;
            onlyLowSpeed=0;
                 
         if(onlyHighSpeed)
            speedFilteringArray=highSpeedIdxes;
         elseif(onlyLowSpeed)
             speedFilteringArray=lowSpeedIdxes;
         end
         
             %}
            allSpeed=1;
            
         allSpeedIdxes= masterSpeedPerTime<maxSpeed & masterSpeedPerTime>minSpeed;
       
         if(allSpeed)
             filteringArray=allSpeedIdxes;
         end
         
         %filteringArray=filteringArray & masterFieldNumPerTime==currFieldNum;
         
         nBins=50;
           %nBins=32;
         %nBins=100;
          inTimeIdxes=~isnan(masterLapTimePerTime) & filteringArray;
         [Ndp,cdp]=hist3([masterDistPerTime(inTimeIdxes),masterSpikePhasePerTime(inTimeIdxes)],[nBins nBins]*1);
         %xlim([0 1])
         masterLapTimePerTime(masterLapTimePerTime>(1-zeroVal))=NaN;
         
         inFieldIdxes=~isnan(masterDistPerTime)& filteringArray;
         
         logTime=log(1-masterLapTimePerTime(inFieldIdxes));
         logTime(logTime<-1.2)=NaN;
         %try
         %[Ntp,ctp]=hist3([masterLapTimePerTime(inFieldIdxes),masterSpikePhasePerTime(inFieldIdxes)],[nBins nBins]*1);
         [Ntp,ctp]=hist3([1-masterLapTimePerTime(inFieldIdxes),masterSpikePhasePerTime(inFieldIdxes)],[nBins nBins]*1);

         %[Ntp,ctp]=hist3([logTime(:),masterSpikePhasePerTime(inFieldIdxes)],[nBins nBins]*1);

         Ntp=flipud(Ntp);
         
         %catch
         %    continue
         %end
         %{
         
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %find log scale circ lin corr of avg phase for each cell per time
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         adj=0;
          %adj=10;
         flippedTime=1-masterLapTimePerTime(inFieldIdxes)+adj;
         timeInFieldPhases=masterSpikePhasePerTime(inFieldIdxes);
         
         figure(logFitFigTime)
         
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %phase vs log time
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           subplot(ceil(sqrt(numSpeedStrata)),ceil(sqrt(numSpeedStrata))+1,speedStratum+1)
           
           xAxisBinEdges=0:0.01:1;
           xAxisBinCenters=edgesToBins(xAxisBinEdges);
           
           avgPhasePerBin=NaN(length(xAxisBinCenters),1);
           for bi=1:length(xAxisBinCenters)
               currBinStart=xAxisBinEdges(bi);
               currBinEnd=xAxisBinEdges(bi+1);
               currBinIdxes=flippedTime>currBinStart & flippedTime<currBinEnd;
               avgPhasePerBin(bi)=circMeanDeg(timeInFieldPhases(currBinIdxes));
               if(circMeanDeg(timeInFieldPhases(currBinIdxes))==0)
                   avgPhasePerBin(bi)=NaN;
               end
           end
        
           sigmoidCenter=0.5;
           sigmoidSteepness=10;
           %sigmoidSteepness=15;
      
           
         plot(flippedTime,timeInFieldPhases,'k.','MarkerSize',15)
          %plot(sigmoid(flippedTime,sigmoidCenter,sigmoidSteepness),timeInFieldPhases,'k.','MarkerSize',15)
         
              %plot(xAxisBinCenters,avgPhasePerBin,'k.','MarkerSize',15)
         maxTime=0.95;
         maxTime=0.8;
          maxTime=0.9;
                  %maxTime=1;
           %maxTime=0.7;
             % maxTime=1;
         timeInFieldPhases(log(flippedTime)<log(1-maxTime))=NaN;
         flippedTime(log(flippedTime)<log(1-maxTime))=NaN; %outliers bias

         [ rhoCircLog,pCircLog,s,b ] = kempter_lincirc( log(flippedTime),timeInFieldPhases/360*2*pi );
         %[ rhoCircLog,pCircLog,s,b ] = kempter_lincirc( sigmoid(flippedTime,sigmoidCenter,sigmoidSteepness),timeInFieldPhases/360*2*pi );

         %rhoCircLog=abs(rhoCircLog);%flips sign randomly???
         %rhoCircLogAllFields=[rhoCircLogAllFields; rhoCircLog];
         rhoCircLogAllFields=[rhoCircLogAllFields; s*360];
         
                %neitherIsNaNIdxes=~isnan(flippedTime) & ~isnan(timeInFieldPhases);
         %[corrCoeffs,ps]=corrcoef([log(flippedTime(neitherIsNaNIdxes)),timeInFieldPhases(neitherIsNaNIdxes)/360]);
          %[corrCoeffs,ps]=corrcoef([(flippedTime(neitherIsNaNIdxes)),timeInFieldPhases(neitherIsNaNIdxes)/360]);
         %R=corrCoeffs(1,2);
         %p=ps(1,2);
         
          %title({sprintf('speed from %.2f to %.2f m/s',minSpeed,maxSpeed),sprintf('Phase vs log(t), circ-linear R = %.2f,p=%.3f',rhoCirc,pCirc)})
          title({'Phase vs log(time to end of field)',sprintf('circ-linear R = %.2f, slope=%.1f',rhoCircLog,s*360)})
           %title({'Phase vs sigmoid(time to end of field)',sprintf('circ-linear R = %.2f',rhoCircLog)})
            xlim([(1-maxTime) 1+adj])
              ylim([0 360])
         %xlim([-1 0])
         set(gca,'xscale','log')
          
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %phase vs linear time
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            subplot(ceil(sqrt(numSpeedStrata)),ceil(sqrt(numSpeedStrata))+1,speedStratum)
                plot(flippedTime,timeInFieldPhases,'k.','MarkerSize',15)
                  %plot(xAxisBinCenters,avgPhasePerBin,'k.','MarkerSize',15)
            
                  numer=sum(flippedTime>0.5 & timeInFieldPhases/360>0.5);
                  denom=sum(flippedTime<0.5 & timeInFieldPhases/360<0.5);
                    asymIdx=numer/denom;
                    allAsymIdx=[allAsymIdx;asymIdx];
            [ rhoCircLin,pCircLin,sLin,b ] = kempter_lincirc( (flippedTime),timeInFieldPhases/360*2*pi );
        
         %rhoCircLin=abs(rhoCircLin);%flips sign randomly???
         %rhoCircLinAllFields=[rhoCircLinAllFields; rhoCircLin];
          rhoCircLinAllFields=[rhoCircLinAllFields; sLin*360];
         
            %[corrCoeffs,ps]=corrcoef([log(flippedTime(neitherIsNaNIdxes)),timeInFieldPhases(neitherIsNaNIdxes)/360]);
          %[corrCoeffs,ps]=corrcoef([(flippedTime(neitherIsNaNIdxes)),timeInFieldPhases(neitherIsNaNIdxes)/360]);
        % R=corrCoeffs(1,2);
        % p=ps(1,2);
         
          %title({sprintf('speed from %.2f to %.2f m/s',minSpeed,maxSpeed),sprintf('Phase vs log(t), circ-linear R = %.2f,p=%.3f',rhoCirc,pCirc)})
          title({'Phase vs time to end of field,',sprintf('circ-linear R = %.2f, slope=%.3f',rhoCircLin,sLin*360)})
          %title({'Phase vs (time to end of field)',sprintf('circ-linear R = %.2f',rhoCircLin)})
          
           %title({sprintf('speed from %.2f to %.2f m/s',minSpeed,maxSpeed),sprintf('R = %.2f,p=%.3f',R,p)})
         xlim([(1-maxTime) 1+adj])
         ylim([0 360])
         %xlim([-1 0])
         %set(gca,'xscale','log')
         maxFig
         setFigFontTo(20)
         drawnow
         %pause(1)
         continue
         %{
     
          figure(logFitFigDist)
         flippedDist=1-masterDistPerTime(inTimeIdxes);
         distInFieldPhases=masterSpikePhasePerTime(inTimeIdxes);
         subplot(ceil(sqrt(numSpeedStrata)),ceil(sqrt(numSpeedStrata))+1,speedStratum)
         plot(flippedDist,distInFieldPhases,'k.','MarkerSize',15)
         
           [ rhoCirc,pCirc,s,b ] = kempter_lincirc( log(flippedDist),distInFieldPhases/360*2*pi );
           
         
         %[corrCoeffs,ps]=corrcoef([log(flippedTime(neitherIsNaNIdxes)),timeInFieldPhases(neitherIsNaNIdxes)/360]);
          %[corrCoeffs,ps]=corrcoef([(flippedTime(neitherIsNaNIdxes)),timeInFieldPhases(neitherIsNaNIdxes)/360]);
         %R=corrCoeffs(1,2);
         %p=ps(1,2);
         
          title({sprintf('speed from %.2f to %.2f m/s',minSpeed,maxSpeed),sprintf('Phase vs log(x), circ-linear R = %.2f,p=%.3f',rhoCirc,pCirc)})
           %title({sprintf('speed from %.2f to %.2f m/s',minSpeed,maxSpeed),sprintf('R = %.2f,p=%.3f',R,p)})
         xlim([(1-maxTime) 1])
         %}
             continue
               end
         
          %%
               figure; 
              %{
            histogram(rhoCircLinAllFields,linspace(-0.5,1,25))
               hold on
                histogram(rhoCircLogAllFields,linspace(-0.5,1,25))
                %}
                 histogram(rhoCircLinAllFields,linspace(-500,500,25))
               hold on
                histogram(rhoCircLogAllFields,linspace(-500,500,25))
          
                [Nall,edges]=histcounts(allAsymIdx,linspace(0,5,20))
          
               
     
                
                bootstrapSampSize=100;
                
                numShuffles=10000;
                
                nonNanTimes=masterLapTimePerTime(~isnan(masterLapTimePerTime));
                nonNanPhases=masterSpikePhasePerTime(~isnan(masterLapTimePerTime));
                
                for si=1:numShuffles
                    si
                    subSampleTimeIdxes=randsample(length(nonNanTimes),bootstrapSampSize);
                    sampleTimes=1-nonNanTimes(subSampleTimeIdxes);
                    
                    %{
                    subSamplePhaseIdxes=randsample(length(nonNanPhases),bootstrapSampSize);
                    samplePhases=nonNanPhases(subSamplePhaseIdxes)/360;
                    %}
                    
                    sampleDegrees=randsample((0:0.1:360),bootstrapSampSize);
                    samplePhases=sampleDegrees/360;
                    
                    
                    sampleNumer=sum(sum(sampleTimes>0.5 & samplePhases>0.5));
                    sampleDenom=sum(sum(sampleTimes<0.5 & samplePhases<0.5));
                    bootstrapAsymIdxes(si)=sampleNumer/sampleDenom;
                end
                hold on
                [Nboot,edgesBoot]=histcounts(bootstrapAsymIdxes,linspace(0,5,50))
                      figure;
                [p1]=plot(edgesToBins(edges),Nall/sum(Nall),'r-','LineWidth',5)
                hold on
                 p2=plot(edgesToBins(edgesBoot),Nboot/sum(Nboot),'k-','LineWidth',5)
                  plot([ 1 1],ylim,'k--','LineWidth',5)
                             xlim([0 5])
                             legend([p2 p1],{'uniform phase, shuffled times',sprintf('%d place fields from %d cells',numFieldsAnalyzed,cellCount)})
                             xlabel('Asymmetry index')
                              ylabel('Probability')
                              %title('Quantifying linearity of real place cell phase precession')
                              maxFig
                              setFigFontTo(24)
                              saveas(gcf,'asymIdx.tif')
                              %}
                             
         
         
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %convert to prob
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
         Ndp=Ndp/sum(Ndp(:));
         Ntp=Ntp/sum(Ntp(:));
         
         [distCount,distEdges]=histcounts(masterDistPerTime,nBins);
         [timeCount,timeEdges]=histcounts(masterLapTimePerTime,nBins);
         [spikeCount,spikeEdges]=histcounts(masterSpikePhasePerTime,nBins);
         
          %time outlier
          %timeCount(1)=NaN;
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %smooth prob
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         movWindSize=10;
         distProb=smooth(distCount/nansum(distCount),movWindSize);
         timeProb=smooth(timeCount/nansum(timeCount),movWindSize); %moving average smooth`
         
         spikeProb=smooth(spikeCount/nansum(spikeCount),movWindSize);
         
         distBins=edgesToBins(distEdges);
         timeBins=edgesToBins(timeEdges);
         
         %{
         figure
         subplot(2,1,1)
         plot(timeBins,timeCount)
         hold on
          plot(timeBins,smooth(timeCount,10))
          subplot(2,1,2)
         plot(distBins,distCount)
         hold on
          plot(distBins,smooth(distCount,10))
         %}
         
        
         %normalize by probablity of column
         for yBin=1:nBins
             Ndp(:,yBin)=Ndp(:,yBin)./distProb(:);
             Ntp(:,yBin)=Ntp(:,yBin)./timeProb(:);
             %{
             if(yBin==nBins)
                 Ntp(:,yBin)=NaN(size(Ntp(:,yBin)));
                 Ndp(:,yBin)=NaN(size(Ndp(:,yBin)));
             end
             %}
                 
         end
         
         for xBins=1:nBins
          %   Ndp(xBins,:)=Ndp(xBins,:)./spikeProb(:)';
           %  Ntp(xBins,:)=Ntp(xBins,:)./spikeProb(:)';
         end
         
         %Ndp=imgaussfilt(Ndp,[.1 .1]*5, 'Padding','circular');
         %Ntp=imgaussfilt(Ntp,[.1 .1]*5, 'Padding','circular');
         %figure; histogram(masterSpikePhasePerTime,120)
         Ndp(end,:)=NaN;
          Ntp(end,:)=NaN;
          
           Ndp(1,:)=NaN;
          Ntp(1,:)=NaN;
          
         Ndp=nanGaussSmoothPhaseProb(Ndp)
         Ntp=nanGaussSmoothPhaseProb(Ntp)
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %divide by prior probability of rat being in that bin
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
         Ndp=Ndp/nansum(Ndp(:));
         Ntp=Ntp/nansum(Ntp(:));
         
         figure(distFig)
         subplot(ceil(sqrt(numSpeedStrata)),ceil(sqrt(numSpeedStrata)),speedStratum)
         %imagesc(flipud(N'))
           Ndp=nanGaussSmoothPhaseProb(Ndp);
         omarPcolor(cdp{1},cdp{2},Ndp',distFig)
          cb=colorbar
         ylabel(cb,'Joint probability')
         colormap(gca,jet)
         xlabel('Dist in field (frac)')
         ylabel('Spike phase (deg)')
           xlim([0 1]-zeroVal)
           pMax=0.08*50/nBins;
           title(sprintf('speed from %.2f to %.2f m/s',minSpeed,maxSpeed))
            caxis([0 prctile(Ndp(:),99)])
            daspect([1 360 1])
           uberTitle('Distance in field phase precession')
              
        setFigFontTo(14)
    maxFig
    saveas(gcf,'DistanceInFieldPhasePrecessionAvg.tif')
            %caxis([0 pMax])
         %figure; 
         %plot(masterDistPerTime,masterSpikePhasePerTime,'k.')
         figure(timeFig)
         subplot(ceil(sqrt(numSpeedStrata)),ceil(sqrt(numSpeedStrata)),speedStratum)
         
         %xlim([0 1])
         
         %imagesc(flipud(N'))
         Ntp=nanGaussSmoothPhaseProb(Ntp);
         omarPcolor(ctp{1},ctp{2},Ntp',timeFig)
         cb2=colorbar
         ylabel(cb2,'Joint probability')
         colormap(gca,jet)
         xlabel('Time in field (frac)')
              %xlabel('log(Time in field)')
         ylabel('Spike phase (deg)')
          %xlim([0 1]-zeroVal)
           title(sprintf('speed from %.2f to %.2f m/s',minSpeed,maxSpeed))
           daspect([1 720 1])
           %axis square
         %caxis([0 prctile(Ntp(:),95)])
           caxis([0 prctile(Ntp(:),99)])
           ylim([50 360])
           xlim([0 0.8])
          uberTitle('Time in field phase precession')
       
        setFigFontTo(14)
    maxFig
    saveas(gcf,'TimeInFieldPhasePrecessionAvg.tif')
    
         end
        %{
         
         dBins=latestOutputDataStruct.distBinsPhase;
         tBins=latestOutputDataStruct.timeBinsPhase;
         p=latestOutputDataStruct.meanPhaseSmooth;
         combinedDTphase=NaN(size(dBins));
         for d=1:length(dBins)
             combinedDTphase(d)=p(d,d);
         end
         %figure
         %plot(combinedDTphase,'ko')
    %}
    
    
end
clef()
%autoArrangeFigures
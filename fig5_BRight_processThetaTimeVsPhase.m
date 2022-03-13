close all; clear all; clc
dataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles';
%xxxxx
setTightSubplots_SpaceTime

filePaths=getFilePathsRegex(dataDir,'*mat');

maxNumSpikesInCycle=500;
maxNumCycles=50;

minNumSpikes=30;
maxCycDisp=20;
maxCycDisp=15;
minSpikeSpeed=0.05;
%minNumPts=150;
minNumPts=100;
minNumPts=50;
%minSpikeSpeed=0.1;

colorCodeForSpeed=0;
colorCodeForPhase=1;
showPlots=0;

maxCorr=-0.05;
maxCorr=-0.1;
maxCorr=-0.2;
%maxCorr=-0.15;

justFirstSpikes=1;
%justFirstSpikes=0;

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
figure
plottedUnitNum=1;
for fi=startFile:length(filePaths)
       fi
       
       disp('loading data')
       tic
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
    toc

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
    
        if(~isfield(data.directionSpecificStats,currFieldDirection))
            continue
        end
        numFields=size(data.directionSpecificStats.(currFieldDirection).expectedWidthPerField,1);

        
        for fii=1:numFields
            disp(sprintf('processing cell %d, direction %d, field%d',fi, di, fii))
           
            tic
            if(~isfield(data.thetaBasedTimeVsPhaseInfo.(currFieldDirection),sprintf('field%d',fii)))
                %currFieldThetaTimeData=[];
                continue
            end
            currFieldThetaTimeData=data.thetaBasedTimeVsPhaseInfo.(currFieldDirection).(sprintf('field%d',fii));
            
            expectedNumThetaCyclesThisField=median(currFieldThetaTimeData.numThetaCyclesInFieldPerLap);
            
            spikeCyclesFromFieldStart=currFieldThetaTimeData.allLapSpikeCycleIDsFromFieldEntry;
        
            spikePhases=currFieldThetaTimeData.allLapInFieldSpikePhases;
            spikeSpeeds=currFieldThetaTimeData.allLapInFieldSpikeSpeeds;
            
            if(sum(spikeCyclesFromFieldStart<maxCycDisp)<minNumSpikes)
                currFieldThetaTimeData=[];
                continue
            end
            
             if(fi==7)
                disp('')
                continue %an obvious repeat of fi=6...
             end
             
              if(fi==24)
                disp('')
                continue %an obvious repeat of fi=23...
              end
              
              if(fi==73)
                disp('')
                continue %an obvious repeat of earlier...
              end
             
            nanIdxes=isnan(spikeCyclesFromFieldStart) | isnan(spikePhases);
            
            domainIdxes=spikeCyclesFromFieldStart<maxCycDisp;
            
            goodSpeedIdxes=spikeSpeeds>minSpikeSpeed;
            
            goodIdxes=domainIdxes & ~nanIdxes & goodSpeedIdxes;
            numPts=sum(goodIdxes);
            spikeCyclesFromFieldStart=spikeCyclesFromFieldStart(goodIdxes);
            spikePhases=spikePhases(goodIdxes);
            spikeSpeeds=spikeSpeeds(goodIdxes);
            
             [m,b,R]=getLinearFit(spikeCyclesFromFieldStart,spikePhases);
             %[ rho,p,s,b ]=kempter_lincirc(spikeCyclesFromFieldStart(~nanIdxes),ang2rad( spikePhases(~nanIdxes)));
            if(R>=maxCorr)
                continue
            end
           
            if(numPts<minNumPts)
                continue
            end
            
            if(min(spikeCyclesFromFieldStart)>4)
                continue
                
            end
              

            %spikePhasesPerCycleNum=NaN(maxNumSpikesInCycle,maxNumCycles);
            
            circMeanSpikePhasePerCycleNum=NaN(maxNumCycles,1);
            circKappaSpikePhasePerCycleNum=NaN(maxNumCycles,1);
            
            for ci=1:maxNumCycles
                spikePhasesInThisCycle=spikePhases(spikeCyclesFromFieldStart==ci);
                %numPhasesInThisCycle=length(spikePhasesInThisCycle);
                if(isempty(spikePhasesInThisCycle))
                    continue
                end
                circMeanSpikePhasePerCycleNum(ci)=circMeanDeg(spikePhasesInThisCycle);
                circKappaSpikePhasePerCycleNum(ci)=circ_kappa(ang2rad(spikePhasesInThisCycle));
                %spikePhasesPerCycleNum(1:numPhasesInThisCycle,currCycleNum)=spikePhasesInThisCycle;
            end 
            toc
            
            %subplot(6,9,plottedUnitNum)
             %subplot(7,7,plottedUnitNum)
             %subplot(7,8,plottedUnitNum)
             %subplot(6,8,plottedUnitNum)
             if(fi>=10)
             figure
            %plot(fractionOfFieldTime,circMeanSpikePhasePerCycleNum)
            
            dispHighSpeed=0.6;
            dispLowSpeed=0.3;
            
            dispHighSpeed= prctile(spikeSpeeds,80);
            dispLowSpeed=prctile(spikeSpeeds,20);
           
            highSpeedSpikeIdxes=spikeSpeeds>=dispHighSpeed;
             lowSpeedSpikeIdxes=spikeSpeeds<=dispLowSpeed;
            
            %plot(spikeCyclesFromFieldStart,spikePhases,'k.','MarkerSize',5)
            markerSize=35;
            plot(spikeCyclesFromFieldStart(highSpeedSpikeIdxes),spikePhases(highSpeedSpikeIdxes),'r.','MarkerSize',markerSize)
            
            box off
            hold on
            
            plot(spikeCyclesFromFieldStart(highSpeedSpikeIdxes),spikePhases(highSpeedSpikeIdxes)+360,'r.','MarkerSize',markerSize)
                plot(spikeCyclesFromFieldStart(lowSpeedSpikeIdxes),spikePhases(lowSpeedSpikeIdxes)+360,'b.','MarkerSize',markerSize)
            plot(spikeCyclesFromFieldStart(lowSpeedSpikeIdxes),spikePhases(lowSpeedSpikeIdxes),'b.','MarkerSize',markerSize)
       
                plot(spikeCyclesFromFieldStart(highSpeedSpikeIdxes),spikePhases(highSpeedSpikeIdxes)-360,'r.','MarkerSize',markerSize)
                plot(spikeCyclesFromFieldStart(lowSpeedSpikeIdxes),spikePhases(lowSpeedSpikeIdxes)-360,'b.','MarkerSize',markerSize)
             
           
            %{
                            plot(spikeCyclesFromFieldStart(highSpeedSpikeIdxes),spikePhases(highSpeedSpikeIdxes)+360,'ro','MarkerSize',10)
                plot(spikeCyclesFromFieldStart(lowSpeedSpikeIdxes),spikePhases(lowSpeedSpikeIdxes)+360,'bo','MarkerSize',10)
            plot(spikeCyclesFromFieldStart(lowSpeedSpikeIdxes),spikePhases(lowSpeedSpikeIdxes),'bo','MarkerSize',10)
       
                plot(spikeCyclesFromFieldStart(highSpeedSpikeIdxes),spikePhases(highSpeedSpikeIdxes)-360,'ro','MarkerSize',10)
                plot(spikeCyclesFromFieldStart(lowSpeedSpikeIdxes),spikePhases(lowSpeedSpikeIdxes)-360,'bo','MarkerSize',10)
            %}
            %plot(1:maxCycDisp,circMeanSpikePhasePerCycleNum(1:maxCycDisp),'r-')
            plottedUnitNum=plottedUnitNum+1;
            %axis square
            ylim([0 360])
            ylim([-360 720])
              estMaxNumCycles=prctile(currFieldThetaTimeData.numThetaCyclesInFieldPerLap,85);
              
        %xlim([0 min(30,estMaxNumCycles)])
        %xlim([0 min(maxCycDisp,estMaxNumCycles)])
         xlim([0 maxCycDisp])
         
         daspect([15/1.5 360 1])
         if((plottedUnitNum-1)>=41)
          xlabel('Cycle no. in field')
         end
          if(mod((plottedUnitNum-1),8)==1)
        ylabel('Phase (deg)')
          end
          
            xlabel('Cycle no. in field')
                ylabel('Phase (deg)')
                drawnow
                disp('')
            end

            %title(sprintf('R=%.2f',R))
            %title(sprintf('num pts=%.2f',numPts))
            %meanKappa=NaN;
        end
        %meanKappa=nanmean(circKappaSpikePhasePerCycleNum);
      
        
         
            
            
        %{
        if(nanmean(circKappaSpikePhasePerCycleNum)>1)
            continue
        end
        %}
        %circMeanSpikePhasePerCycleNum(circKappaSpikePhasePerCycleNum>2)=NaN;
        
        %fractionOfFieldTime=(1:maxNumCycles)/expectedNumThetaCyclesThisField;
        
        %{
        if(mean(circMeanSpikePhasePerCycleNum(1:floor(expectedNumThetaCyclesThisField/2)))<180)
            continue
        end
        %}
        %{
        subplot(10,10,plottedUnitNum)
        %plot(fractionOfFieldTime,circMeanSpikePhasePerCycleNum)
        plot(circMeanSpikePhasePerCycleNum)
        hold on
        %}
      
        %try
            %xlim([0 min(50,expectedNumThetaCyclesThisField*1.5)])
        %end
        
      
        %{
        subplot(2,1,2)
        plot(fractionOfFieldTime,circKappaSpikePhasePerCycleNum)
        hold on
        ylim([0 10])
        %}
        %maxFigHalfWidth
        maxFigMukkaalWidth
        setFigFontTo(12)
        drawnow
        
   end
    
end

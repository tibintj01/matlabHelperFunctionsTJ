function [] = plotTimeFieldAndThetaPhases(currUnitInfo,saveDir)

for di=1:2
    
    if(di==1)
        currDirectionStr='rightward';
    else
         currDirectionStr='leftward';
    end
    
    if(~isfield(currUnitInfo.directionSpecificStats,currDirectionStr))
        continue
    end
  
    if(~isfield(currUnitInfo,sprintf('manualFieldStartsM%s',currDirectionStr)))
        continue
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %load variables from unit info struct
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    positionTimeAxis=currUnitInfo.positionTimeAxis;
    fieldEntryTimePerTime=currUnitInfo.directionSpecificStats.(currDirectionStr).fieldEntryTimePerTime;
    fieldNumPerTime=currUnitInfo.directionSpecificStats.(currDirectionStr).fieldNumPerTime;
     lapNumPerTime=currUnitInfo.directionSpecificStats.(currDirectionStr).lapNumberPerTime;
    
    fieldStarts=currUnitInfo.(sprintf('manualFieldStartsM%s',currDirectionStr));
    
    
    
    numFields=size(currUnitInfo.directionSpecificStats.(currDirectionStr).fieldWidthsPerLap,1);
    
    
    
    if(di==1)
        currSpikePhases=currUnitInfo.rightSpikePhases;
        currSpikeTimes=currUnitInfo.rightSpikeTimes;
    else
        currSpikePhases=currUnitInfo.leftSpikePhases;
        currSpikeTimes=currUnitInfo.leftSpikeTimes;
    end
    
    spikeTimesFromLatestFieldEntry=NaN(numFields,length(currSpikeTimes));
    lapNumPerSpike=NaN(numFields,length(currSpikeTimes));
    numSpikes=length(currSpikeTimes);
    
    if(length(positionTimeAxis)~=length(fieldEntryTimePerTime))
        continue
    end

    for si=1:numSpikes
        currSpikeTime=currSpikeTimes(si);
        [~,closestPrevPosTimeIdx]=min(abs(currSpikeTime-positionTimeAxis));
        if(positionTimeAxis(closestPrevPosTimeIdx)>currSpikeTime)
            closestPrevPosTimeIdx=closestPrevPosTimeIdx-1;
        end
        currSpikeFieldEntryTime=fieldEntryTimePerTime(closestPrevPosTimeIdx);
        currSpikeLapNum=lapNumPerTime(closestPrevPosTimeIdx);

        for fi=1:numFields
            %closestPrevPosTime=positionTimeAxis(closestPrevPosTimeIdx);
           if(fi~= fieldNumPerTime(closestPrevPosTimeIdx))
               continue
           end
            spikeTimesFromLatestFieldEntry(fi, si)=currSpikeTime-currSpikeFieldEntryTime;
            lapNumPerSpike(fi,si)=currSpikeLapNum;
        end
    end
   

    rightCycleAmps=currUnitInfo.cycleAmpPerRightSpike;
    leftCycleAmps=currUnitInfo.cycleAmpPerLeftSpike;    
    
    ampDispThresh=2.5;
    rightCycleAmps(rightCycleAmps<ampDispThresh)=NaN;
    leftCycleAmps(leftCycleAmps<ampDispThresh)=NaN;
    
    rightCycleAmps=(rightCycleAmps-ampDispThresh)*7.5;
    leftCycleAmps=(leftCycleAmps-ampDispThresh)*7.5;
    
    if(di==1)
        currCycleAmps=rightCycleAmps;
    else
        currCycleAmps=leftCycleAmps;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %separate different directions into separate variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    speedPerTimeStep=currUnitInfo.speedPerTimeStep;
    rightwardTimeIdxes=speedPerTimeStep>0.05;
    leftwardTimeIdxes=speedPerTimeStep<-0.05;
    
    rightwardSpeedPerTimeStep=speedPerTimeStep;
    leftwardSpeedPerTimeStep=abs(speedPerTimeStep);
    
    rightwardSpeedPerTimeStep(~rightwardTimeIdxes)=NaN;
    leftwardSpeedPerTimeStep(~leftwardTimeIdxes)=NaN;
    
    
    positionPerTimeStep=currUnitInfo.positionPerTimeStep;
    
    baseName=sprintf('%s_6-12Theta_Unit%dInfo_%s',currUnitInfo.unitInfo.sessionName,currUnitInfo.unitInfo.unitIDnum,currUnitInfo.unitInfo.unitCellType);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %plot time fields
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for fi=1:numFields
        if(fi>length(fieldStarts))
            continue
        end
        
        figH=figure;
              
        subplot(2,4,[1 2 3])
        
        if(~iscell(currUnitInfo.imagePaths))
          if(~contains(currUnitInfo.imagePaths,currDirectionStr))
              continue
          end
        imshow(currUnitInfo.imagePaths)
      else
          if(di==1)
              ii=2;
          else
              ii=1;
          end
           if(~contains(currUnitInfo.imagePaths{ii},currDirectionStr))
              continue
          end
         imshow(currUnitInfo.imagePaths{ii})
       end
        
      [durN,durEdges,durBinLabels]=histcounts(spikeTimesFromLatestFieldEntry(fi,:),100);
        smoothWidth=10;
        smoothedTimeCount=smooth(durN,smoothWidth);
        durBins=edgesToBins(durEdges);
        [~,mostCommonDurBin]=max(smoothedTimeCount);
             
        typicalDur=nanmean(spikeTimesFromLatestFieldEntry(fi,durBinLabels==mostCommonDurBin));
        %approxMaxFieldDuration=typicalDur*3;
        
        approxMaxPrctile=prctile(spikeTimesFromLatestFieldEntry(fi,:),95);
        approxMaxFieldDuration=nanmean([typicalDur*3 approxMaxPrctile]);
        
        if(isnan(approxMaxFieldDuration))
            approxMaxFieldDuration=3;
        end
        subplot(2,4,[4])
        histogram(spikeTimesFromLatestFieldEntry(fi,:),100)
        hold on
        plot(durBins,smoothedTimeCount,'k-','LineWidth',5)
        xlabel('Time in field (sec)')
        ylabel('Spike count')
        title({sprintf('Est. max field duration: %.2f sec',approxMaxFieldDuration)})
        
        if(~isnan(currUnitInfo.(sprintf('manualTimeField%dEnd%sSec',fi,currDirectionStr))))
            approxMaxFieldDuration=currUnitInfo.(sprintf('manualTimeField%dEnd%sSec',fi,currDirectionStr));
        end
        
        pR=plot([approxMaxFieldDuration approxMaxFieldDuration],ylim,'r-','LineWidth',2)
        pP=plot([approxMaxPrctile approxMaxPrctile],ylim,'b--','LineWidth',2)
        pT=plot([typicalDur*3 typicalDur*3],ylim,'k--','LineWidth',2)
        
        legend([pT pR pP],'3x typical time','avg','95th percentile time')
        set(gca,'FontSize',12)
        
        subplot(2,1,2)
        %hold on
        %plot(leftSpikePositions,leftSpikePhases,'k.')
        lapNumsThisField=lapNumPerSpike(fi,:);
        scatter(spikeTimesFromLatestFieldEntry(fi,:),currSpikePhases,currCycleAmps*50,lapNumsThisField(:),'.') %'k'
        colormap(jet)
        cbL=colorbar
        ylabel(cbL,'lap no.')
        %xlim([min(positionBins) max(positionBins)])
        xlim([0 approxMaxFieldDuration])
        ylim([0 360])
        ylabel('Theta phase (deg)')
        if(approxMaxFieldDuration<4)
            plotLinesEvery5cm
        else
            plotLinesEveryTenth
        end
        xlabel('Time in field (sec)')
         %box off
         box on
         setFigFontTo(18)
              
         if(approxMaxFieldDuration<2)
              set(gca,'FontSize',18)
         else
            set(gca,'FontSize',12)
         end
        
%spikeTimesFromLatestFieldEntry

title('Theta phase vs time in field')

	set(gca,'XMinorTick','on')
    

        uberTitle(sprintf('%s %s, field %d (starts at %.2f m)',removeUnderscores(baseName),currDirectionStr,fi,fieldStarts(fi)))
        %maxFigManual2d(3,0.8)
        set(gca,'FontSize',18)
        maxFig
        
        

	print('-r150',fullfile(saveDir,sprintf('%s_TimeFieldAndPhaseVsPos_%s_Field%d',baseName,currDirectionStr,fi)),'-dpng') 
       
         %saveas(gcf,fullfile(saveDir,sprintf('%s_TimeFieldAndPhaseVsPos_%s.tif',baseName,saveDirStr)))
           %disp('')
           close all
           
    end%field loop
           
end %direction loop
           %imshow(fullfile(saveDir,sprintf('%s_TimeFieldAndPhaseVsPos_%s.png',baseName,saveDirStr)))
           %disp('')
           
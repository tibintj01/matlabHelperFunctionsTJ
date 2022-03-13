close all; clear all; clc
dataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles';
cycleDataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/';
showPlots=1;

%Cicero_09172014 has only 1 strongly precessing cell used
%sessionNames={'Achilles_10252013','Achilles_11012013','Buddy_06272013', 'Cicero_09102014','Cicero_09172014','Gatsby_08282013'};
sessionNames={'Achilles_10252013','Achilles_11012013','Buddy_06272013', 'Cicero_09102014','Gatsby_08282013'};

%sessionNames={'Achilles_11012013'};
%sessionNames={'Gatsby_08282013'};

halfLineHeight=0.01; %raster line height
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%loop through each session and all corresponding units
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for si=1:length(sessionNames) 
    currSessName=sessionNames{si};
    clearvars -except si currSessName sessionNames showPlots dataDir cycleDataDir
        %{
    currSessMasterUnitData={};
    for fi=1:length(currSessUnitDataPaths)
        currData=load(currSessUnitDataPaths{fi});
        currSessMasterUnitData{fi}=currData;
    end
        %}
    
        if(exist(sprintf('%sSpikeData.mat',currSessName),'file'))
            load(sprintf('%sSpikeData.mat',currSessName))
        else
    [minTime, maxTime,numLapsPerDir]=getSessionTimeBounds(sessionNames{si});
    %maxTime=minTime+10;
    
    [leftwardRefUnitData,rightwardRefUnitData]=getRefUnitData(currSessName);
    
    currSessUnitDataPaths=getRegexFilePaths(dataDir,sprintf('%s*Unit*.mat',currSessName));
        numUnitsInSess=length(currSessUnitDataPaths);
            currSessMasterAllSpikeData={};
            currSessMasterFieldMiddleData={};
            currSessMasterEnterFieldTimeData={};
            %currSessMasterFirstSpikeData={};

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%for each direction, loop through all units
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
                            allTimesStr=(sprintf('%sSpikeTimes',currDir));
                            firstTimesStr=(sprintf('%sSpikeTimesOnlyRefCycleFirsts',currDir));
                            
                            if(isfield(currUnitData.directionSpecificStats,currDirStr))
                            %fieldMiddlePosPerTime=(currUnitData.directionSpecificStats.(currDirStr).fieldExitPosPerTime + ...
                             %   currUnitData.directionSpecificStats.(currDirStr).fieldEntryPosPerTime)/2;
                                fieldNumPerTime=currUnitData.directionSpecificStats.(currDirStr).fieldNumPerTime;
                                fieldEntryTimePerTime=currUnitData.directionSpecificStats.(currDirStr).fieldEntryTimePerTime;
                                fieldExitTimePerTime=currUnitData.directionSpecificStats.(currDirStr).fieldExitTimePerTime;
                                
                                numFields=max(fieldNumPerTime);
                                fieldCenters=(currUnitData.(sprintf('manualFieldStartsM%s',currDirStr))+currUnitData.((sprintf('manualFieldEndsM%s',currDirStr))))/2;
                                 fieldMiddlePosPerTime=NaN(size(fieldNumPerTime));
                                for ti=1:length(fieldNumPerTime)
                                    if(~isnan(fieldNumPerTime(ti)))
                                     fieldMiddlePosPerTime(ti)=fieldCenters(fieldNumPerTime(ti));
                                    end
                                end
                            else
                                currSessMasterAllSpikeData{fi,di}=NaN;
                                currSessMasterFieldMiddleData{fi,di}=NaN;
                                currSessMasterEnterFieldTimeData{fi,di}=NaN;
                                currSessMasterExitFieldTimeData{fi,di}=NaN;
                                continue
                                
                            end
                            if(isfield(currUnitData,allTimesStr))
                                currSessMasterAllSpikeData{fi,di}=currUnitData.(allTimesStr);
                            else
                                currSessMasterAllSpikeData{fi,di}=NaN;
                            end
                            
                            currSessMasterFieldMiddleData{fi,di}=fieldMiddlePosPerTime;
                            currSessMasterEnterFieldTimeData{fi,di}=fieldEntryTimePerTime;
                            currSessMasterExitFieldTimeData{fi,di}=fieldExitTimePerTime;
                            %{
                            if(isfield(currUnitData,firstTimesStr))
                                currSessMasterFirstSpikeData{fi,di}=currUnitData.(firstTimesStr);
                            else
                                currSessMasterFirstSpikeData{fi,di}=NaN;
                            end
                            %}
                end
                toc
            end
            disp('')
            save(sprintf('%sSpikeData',currSessName))
        end


refUnitData=leftwardRefUnitData;
    if(~isstruct(refUnitData))
            continue
        end

    refChCycleInfo=load(refUnitData.refThetaChCyclePath);

    cycleIsDuringMaze=refChCycleInfo.cycleMinTimes>=minTime & refChCycleInfo.cycleMinTimes<maxTime;

    cycleStartTimes=refChCycleInfo.cycleMaxTimes(1:(end-1));
    cycleEndTimes=refChCycleInfo.cycleMaxTimes(2:end);

    cycleInMazeStartTimes=cycleStartTimes(cycleIsDuringMaze);
    cycleInMazeEndTimes=cycleEndTimes(cycleIsDuringMaze);

    numCycles=length(cycleInMazeStartTimes);

    disp('')
    
    filledSignedSpeedPerTime=refUnitData.filledSignedSpeedPerTime;
    positionTimeAxis=refUnitData.positionTimeAxis;
    minSpeed=0.05;
    
     allSpikeTimesByCycleAndUnit={};
     avgSpikeTimesByCycleAndUnit=NaN(numCycles,numUnitsInSess,2);
     unitCenterPerCycleAndField=NaN(numCycles,numUnitsInSess,2);
     unitEntryTimePerCycleAndField=NaN(numCycles,numUnitsInSess,2);
     unitExitTimePerCycleAndField=NaN(numCycles,numUnitsInSess,2);
     speedPerCycle=NaN(numCycles,1);
     timeOfCycleTrough=NaN(numCycles,1);
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %loop through each theta cycle 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic
 
    for ci=1:numCycles

        currCycleStartTime=cycleInMazeStartTimes(ci);
        currCycleEndTime=cycleInMazeEndTimes(ci);
        cycleCenterTime=(currCycleStartTime+currCycleEndTime)/2;
        
        currSpeed=interp1(positionTimeAxis,filledSignedSpeedPerTime,cycleCenterTime);
        speedPerCycle(ci)=abs(currSpeed);
        timeOfCycleTrough(ci)=(currCycleStartTime+currCycleEndTime)/2;
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %SPEED THRESHOLD
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if(currSpeed>minSpeed)
            di=1; %rightward movement
        elseif(currSpeed<-minSpeed)
            di=2;%leftward movement
        else
            allSpikeTimesByCycleAndUnit{ci,fi}=NaN;
            continue
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %find direction during this cycle
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %for di=1:2
        if(di==1)
            currDirStr='rightward';
            currDir='right';
        else
            currDirStr='leftward';
            currDir='left';
        end

    positionTimeAxis=refUnitData.positionTimeAxis;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %loop through each unit in this session, gather spikes in
        %current cycle
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for fi=1:numUnitsInSess
            currUnitDirAllSpikeTimes=currSessMasterAllSpikeData{fi,di};
            currUnitDirAllFieldMiddles=currSessMasterFieldMiddleData{fi,di};
            
            currUnitDirAllFieldEntryTimes=currSessMasterEnterFieldTimeData{fi,di};
            currUnitDirAllFieldExitTimes=currSessMasterExitFieldTimeData{fi,di};
            
            if(isnan(currUnitDirAllFieldMiddles) )
                continue
            end
            
            if(isnan(currUnitDirAllSpikeTimes))
                continue
            end
            
            
            inCycleTimeIdxes=positionTimeAxis>=currCycleStartTime & positionTimeAxis<currCycleEndTime;
            currUnitDirCycleFieldMids=currUnitDirAllFieldMiddles(inCycleTimeIdxes);
            currUnitDirCycleFieldEntryTimes=currUnitDirAllFieldEntryTimes(inCycleTimeIdxes);
            currUnitDirCycleFieldExitTimes=currUnitDirAllFieldExitTimes(inCycleTimeIdxes);
            
            currUnitDirCycleSpikeTimes=currUnitDirAllSpikeTimes(currUnitDirAllSpikeTimes>=currCycleStartTime & currUnitDirAllSpikeTimes<currCycleEndTime);

            %allSpikeTimesByCycleAndUnit{ci,fi}=currUnitDirCycleSpikeTimes-currCycleStartTime;
            avgSpikeTimesByCycleAndUnit(ci,fi,di)=nanmean(currUnitDirCycleSpikeTimes-currCycleStartTime);
                        %avgSpikeTimesByCycleAndUnit(ci,fi,di)=nanmean(min(currUnitDirCycleSpikeTimes))-currCycleStartTime;

            unitCenterPerCycleAndField(ci,fi,di)=nanmean(currUnitDirCycleFieldMids);
            
            unitEntryTimePerCycleAndField(ci,fi,di)=nanmean(currUnitDirCycleFieldEntryTimes);
            unitExitTimePerCycleAndField(ci,fi,di)=nanmean(currUnitDirCycleFieldExitTimes);
        end
        
        if(mod(ci,1000)==0)
            disp(ci/numCycles)
        end
    end
    
   
    toc
    
    minNumCellsPerSeq=5;
    minNumCellsPerSeq=1;
    onlyInFieldSpiking=1;
    
    if(onlyInFieldSpiking)
        avgSpikeTimesByCycleAndUnit(isnan(unitCenterPerCycleAndField))=NaN;
    end
    
    
    if(showPlots)
        fH=figure(11);
        fH2=figure(12);

        sortedMeanCenterIdxPerDi=NaN(size(unitCenterPerCycleAndField,2),2);
        validCycleIdxesPerDi=false(size(unitCenterPerCycleAndField,1),2);
        for dii=1:2
             meanUnitCenters=nanmean(unitCenterPerCycleAndField(:,:,dii),1);
        [~,sortedMeanCenterIdx]=sort(meanUnitCenters);
            sortedMeanCenterIdxPerDi(:,dii)=sortedMeanCenterIdx(:);
             

            for tii=1:size(unitCenterPerCycleAndField,1)
                %if(~isnan(max(unitCenterPerCycleAndField(tii,:,dii))))

                if(sum(~isnan((avgSpikeTimesByCycleAndUnit(tii,:,dii))))>minNumCellsPerSeq)
                    validCycleIdxesPerDi(tii,dii)=true;
                end
            end

            numValidCycles=sum(validCycleIdxesPerDi(:,dii));
            if(numValidCycles==0)
                continue
            end
          
                
            if(dii==1)
                dispH=fH;
            else
                dispH=fH2;
            end
             figure(dispH)
            ax1=subplot(2,1,1)
            omarPcolor(1:numValidCycles,1:numUnitsInSess,squeeze(avgSpikeTimesByCycleAndUnit(validCycleIdxesPerDi(:,dii),sortedMeanCenterIdxPerDi(:,dii),dii))',dispH)
            colormap(jet)
            cb=colorbar
            caxis([0 0.12])
            axis tight
                    xlabel('Theta cycle no.')
            ylabel('Place cell no.')
            ylabel(cb,'Avg spike time in theta cycle (sec)')
            title('Theta sequences in each cycle')

            ax2=subplot(2,1,2)
            omarPcolor(1:numValidCycles,1:numUnitsInSess,squeeze(unitCenterPerCycleAndField(validCycleIdxesPerDi(:,dii),sortedMeanCenterIdxPerDi(:,dii),dii))',dispH)
            colormap(jet)
            cb2=colorbar
            axis tight
            xlabel('Theta cycle no.')
            ylabel('Place cell no.')
            ylabel(cb2,'Field center (m)')
            linkaxes([ax1 ax2],'xy')
            title('Place fields')

            maxUnit=sum(~isnan(meanUnitCenters))+0.5;
            if(~isnan(maxUnit))
                ylim([0.5 maxUnit])
            end
            uberTitle(removeUnderscores(currSessName))
            setFigFontTo(18)
            maxFig
            %close all

        end
    end
   
    save(sprintf('perUnitPerCycleTimingData%s.mat',currSessName),'unitCenterPerCycleAndField','unitEntryTimePerCycleAndField','unitExitTimePerCycleAndField','avgSpikeTimesByCycleAndUnit','numCycles','onlyInFieldSpiking','speedPerCycle','sortedMeanCenterIdxPerDi','validCycleIdxesPerDi','timeOfCycleTrough')
end

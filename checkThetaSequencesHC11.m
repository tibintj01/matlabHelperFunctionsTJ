
close all; clear all; clc
dataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles';
cycleDataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/';
showPlots=1;


%Cicero_09172014 has only 1 strongly precessing cell used
sessionNames={'Achilles_10252013','Achilles_11012013','Buddy_06272013', 'Cicero_09102014','Cicero_09172014','Gatsby_08282013'};
%sessionNames={'Achilles_10252013'};

%sessionNames={'Achilles_11012013'};
%sessionNames={'Gatsby_08282013'};

halfLineHeight=0.01; %raster line height

startSi=1;
%startSi=6;
for si=startSi:length(sessionNames)
    
    currSesCycleFilePaths=getRegexFilePaths(cycleDataDir,sprintf('%s_Ch*_6-12_minmax.mat',sessionNames{si}));
    
    currSessName=sessionNames{si};
    [minTime, maxTime,numLapsPerDir]=getSessionTimeBounds(sessionNames{si});
    maxTime=minTime+10;
    
    trackFilePath=sprintf('%sTrackProps.mat',currSessName);
    currTrackData=load(trackFilePath);
    positionTimeAxis=currTrackData.positionTimeAxis;
    
    [leftwardRefUnitData,rightwardRefUnitData]=getRefUnitData(currSessName);
    
    for di=1:2
        %{
        numLaps=numLapsPerDir(di);
        
        if(isnan(numLaps))
            continue
        end
        %}
        figure
        if(di==1)
            currDirStr='rightward';
            refUnitData=rightwardRefUnitData;
            
        else
            currDirStr='leftward';
            refUnitData=leftwardRefUnitData;
        end
        if(~isstruct(refUnitData))
            continue
        end
        lapStartTimes=refUnitData.lapStartTimesPerDir.(currDirStr);
        lapEndTimes=refUnitData.lapStopTimesPerDir.(currDirStr);
        numLaps=length(lapStartTimes)
        
        typicalLapDur=median(lapEndTimes-lapStartTimes)*1.2;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %place sorted fields
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if(~isfield(currTrackData.trackPopInfo,currDirStr))
            continue
        end
        allCurrFieldsStruct=currTrackData.trackPopInfo.(currDirStr);
        allCurrFieldNames=fieldnames(allCurrFieldsStruct);
        numFields=length(allCurrFieldNames);
        
        
        fieldPosCenters=[];
        fieldPosStarts=[];
        fieldPosEnds=[];
        
        tickColors=jet(numFields);
        for fii=1:numFields
            fieldPosCenters=[fieldPosCenters allCurrFieldsStruct.(allCurrFieldNames{fii}).fieldPosCenterM];
            fieldPosStarts=[fieldPosStarts allCurrFieldsStruct.(allCurrFieldNames{fii}).fieldPosStart];
            fieldPosEnds=[fieldPosEnds allCurrFieldsStruct.(allCurrFieldNames{fii}).fieldPosEnd];
        end
        
        minFieldStart=min([fieldPosEnds(:); fieldPosStarts(:)]);
        maxFieldEnd=max([fieldPosEnds(:); fieldPosStarts(:)]);
        
        [sortedFieldCenters, sortedCenterIdx]=sort(fieldPosCenters);
        for li=1:numLaps
            allMinTimes=NaN(numFields,1);
            
            for fii=1:numFields
                %currFieldSpikeLaps=allCurrFieldsStruct.(allCurrFieldNames{sortedCenterIdx(fii)}).inFieldSpikeLapNums;

                currFieldSpikeTimes=allCurrFieldsStruct.(allCurrFieldNames{sortedCenterIdx(fii)}).inFieldSpikeTimes;
                currFieldSpikeLaps=getLapPerSpikeTime(currFieldSpikeTimes,lapStartTimes,lapEndTimes);

                  inFieldAndLapTimes=currFieldSpikeTimes(currFieldSpikeLaps==li);
                  if(~isempty(inFieldAndLapTimes))
                      
                 
                      minTime=min(inFieldAndLapTimes);
                    allMinTimes(fii)=minTime;
                 end
            end
            allFieldsMinTime=min(allMinTimes);
            %allFieldsMinTime=19648-290;
            
            for fii=1:numFields
          
                %currFieldSpikeLaps=allCurrFieldsStruct.(allCurrFieldNames{sortedCenterIdx(fii)}).inFieldSpikeLapNums;
                
                currFieldSpikeTimes=allCurrFieldsStruct.(allCurrFieldNames{sortedCenterIdx(fii)}).inFieldSpikeTimes;
                
                currFieldSpikeLaps=getLapPerSpikeTime(currFieldSpikeTimes,lapStartTimes,lapEndTimes);
               
                
                startTick=fieldPosStarts(sortedCenterIdx(fii));
                endTick=fieldPosEnds(sortedCenterIdx(fii));
                
                if(li>ceil(sqrt(numLaps))*floor(sqrt(numLaps)))
                    continue
                end
                subplot(ceil(sqrt(numLaps)),floor(sqrt(numLaps)),li)
                
                inFieldAndLapTimes=currFieldSpikeTimes(currFieldSpikeLaps==li);
               normTimeFieldAndLapSpikeTimes=inFieldAndLapTimes-allFieldsMinTime;
                if(isempty(inFieldAndLapTimes))
                    continue
                end
                              
                tickColor=tickColors(fii,:);
                %plotRasterStyle(inFieldAndLapTimes-allFieldsMinTime,fieldPosCenters(sortedCenterIdx(fii)),startTick,endTick,tickColor)
                                plotRasterStyle(inFieldAndLapTimes-allFieldsMinTime,fii,NaN,NaN,[0 0 0])

                   hold on
            end
                
                axis tight
                xlabel('Time into lap (sec)')
                %ylabel('Field bounds (m)')
                      ylabel('Field no.')
                title(sprintf('Lap %d',li))
                xlim([0 typicalLapDur*2/3])
                %xlim([0 3])
                  %xlim([10 20])
                  ylim([0 numFields])
                  box off
                %ylim([minFieldStart maxFieldEnd])
                
                   %daspect([1 1 1])
                drawnow
                
         
        end
        maxFig
        %uberTitle(removeUnderscores(sprintf('All precessing cell firing per lap, %d fields, sorted by position, %s, %s', numFields,currSessName,currDirStr)))
        setFigFontTo(14)
     
        saveas(gcf,sprintf('allFieldSequencesRaster%s_%s.tif',currSessName,currDirStr))
         close all
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %LFP cycle based per channel theta sequences
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %for ci=1:length(currSesCycleFilePaths)
       
        %{
        currChannelCycleData=load(currSesCycleFilePaths{ci});
        minPts=currChannelCycleData.cycleMinAmp;
        maxPts=currChannelCycleData.cycleMaxAmp;
        
        minTimes=currChannelCycleData.cycleMinTimes;
        maxTimes=currChannelCycleData.cycleMaxTimes;
        
        allTimes=[minTimes(:); maxTimes(:)];
        allPts=[minPts(:); maxPts(:)];
        %plot(allTimes,allPts,'-')
        
        allTimes=allTimes(allTimes>=minTime & allTimes<=maxTime);
        plotRasterStyle(allTimes,ci)
        
        
         
        hold on
        xlim([minTime, maxTime])
        disp('')
        
        if(mod(ci,8)==0)
            plot(xlim,[ci ci],'k-','LineWidth',5)
        end
        
        
        
    end
    ylim([0 length(currSesCycleFilePaths)+0.5])
    xlabel('Time (sec)')
    ylabel('Ch number')
  
     title(removeUnderscores(sprintf('Theta peaks and trough times, %s',currSessName)))
        setFigFontTo(18)
        
      %}
    
end
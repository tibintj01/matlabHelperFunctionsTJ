close all; clear all; clc
dataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles';
cycleDataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/';
showPlots=1;
zScorePlot=1;
%zScorePlot=0;

minSpeed=0.05;

cycleHalf=0.06; %60 msec
cycleDur=0.12;


plotVert=1;

maxSpaceDiff=1;
maxSpaceDiff=0.5;
%maxSpaceDiff=0.3;
maxSpaceDiff=1;
spaceFracMax=0.3*0.3;
spaceFracMax=0.15;
spaceFracMax=0.5;
spaceFracMax=0.2;
%spaceFracMax=0.15;
%spaceFracMin=0.05;
spaceFracMin=0;
spaceFracMin=0.01;
spaceFracMin=0;
%spaceFracMin=0.05;
%spaceFracMax=0.05;

%{
minSpeedForThetaDiff=0.2;
maxSpeedForThetaDiff=0.6;
%}



minSpeedForThetaDiff=0;
maxSpeedForThetaDiff=10;


lastFrac=0.35;
lastFrac=0.33333;
%lastFrac=0.4;
lastFrac=0.5;
%lastFrac=0.4;

%lastFrac=0.2
cofiringThresh=0.2;
cofiringThresh=0.1;
%cofiringThresh=0.08;

lateCycleTimeThresh=cycleDur*(1-lastFrac);
earlyCycleTimeThresh=cycleDur*lastFrac;
numBins=30;
numBins=50;
%{
lateCycleTimeThresh=cycleDur*0.6;
earlyCycleTimeThresh=cycleDur*0.4;
%}

%{
lateCycleTimeThresh=cycleDur*0.5;
earlyCycleTimeThresh=cycleDur*0.5;
%}

%Cicero_09172014 has only 1 strongly precessing cell used
%sessionNames={'Achilles_10252013','Achilles_11012013','Buddy_06272013', 'Cicero_09102014','Cicero_09172014','Gatsby_08282013'};
sessionNames={'Achilles_10252013','Achilles_11012013','Buddy_06272013', 'Cicero_09102014','Gatsby_08282013'};

%sessionNames={'Achilles_11012013'};
%sessionNames={'Gatsby_08282013'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%loop through each session and load per cycle spike time data across units
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
masterEarlyDiffs=[];
masterLateDiffs=[];
numEarlyPhasePairs=0;
numLatePhasePairs=0;

startSi=1;
%startSi=3;
 allSpeedVsTimeDiffR=[];
for si=startSi:length(sessionNames) 
    currSessName=sessionNames{si};
    currSeqDataFilePath=sprintf('perUnitPerCycleTimingData%s.mat',currSessName);
    
    currThetaSeqData=load(currSeqDataFilePath);
     numUnits=size(currThetaSeqData.avgSpikeTimesByCycleAndUnit,2);
     numCycles=currThetaSeqData.numCycles;
     sortedMeanCenterIdxPerDi=currThetaSeqData.sortedMeanCenterIdxPerDi;
     timeOfCycleTrough=currThetaSeqData.timeOfCycleTrough;
     
     numCyclesCofiring=zeros(numUnits,numUnits,2);
       
numEarlyPhasePairs=0;
numLatePhasePairs=0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %for each direction, loop through each pair of units and each cycle 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for di=1:2
        
        allPlottedLateDiffs=[];
        
        if(exist(sprintf('plottedLateDiffs%s_%d.mat',currSessName,di),'file'))
            sortIdxData=load(sprintf('plottedLateDiffs%s_%d.mat',currSessName,di));
            sortByLateDiffIdx=sortIdxData.sortByLateDiffIdx;
            allPlottedLateDiffsSortedSaved=sortIdxData.allPlottedLateDiffsSorted;
        end
        if(di==1)
            currDiStr='rightward';
            
            if(strContainsCircSessionName(currSessName))
                continue
            end
        else
            currDiStr='leftward';
        end
        timeDifferencesPerCycle=NaN(numUnits,numUnits,numCycles);
        currTimeInMetaFieldPerCycle=NaN(numUnits,numUnits,numCycles);
        numCyclesThreshPerRow=NaN(numUnits,1);
        numCyclesMaxPerRow=NaN(numUnits,1);
        alreadyProcessed=false(numUnits,numUnits);
        
        currFieldCentersPerCycle=currThetaSeqData.unitCenterPerCycleAndField(:,:,di);
        currFieldStartTimesPerCycle=currThetaSeqData.unitEntryTimePerCycleAndField(:,:,di);
        currFieldEndTimesPerCycle=currThetaSeqData.unitExitTimePerCycleAndField(:,:,di);
        
        %loop through all pairs of units
        for unit1Num=1:numUnits
            for unit2Num=1:numUnits %vs upper triangle matrix
            %for unit2Num=(unit1Num+1):numUnits %vs upper triangle matrix
            
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %get timing difference within cycle, to compare with speed,
                %spatial difference, and time in meta-field
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if(alreadyProcessed(unit1Num,unit2Num))
                    continue
                end
              
                for ci=1:numCycles
                    unit1TimingCurrCycle=currThetaSeqData.avgSpikeTimesByCycleAndUnit(ci,unit1Num,di);
                    unit2TimingCurrCycle=currThetaSeqData.avgSpikeTimesByCycleAndUnit(ci,unit2Num,di);
                    
                    %timeFromUnit1Start=timeOfCycleTrough(ci)-currFieldStartTimesPerCycle(ci,unit1Num);
                    %timeFromUnit2Start=timeOfCycleTrough(ci)-currFieldStartTimesPerCycle(ci,unit2Num);
                    
                    timeFromUnit1Start=timeOfCycleTrough(ci)-currFieldStartTimesPerCycle(ci,unit1Num);
                    timeFromUnit2Start=timeOfCycleTrough(ci)-currFieldStartTimesPerCycle(ci,unit2Num);
                    
                    %currTimeInMetaFieldPerCycle(unit1Num,unit2Num,ci)=max(timeFromUnit1Start,timeFromUnit2Start);
                    currTimeInMetaFieldPerCycle(unit1Num,unit2Num,ci)=timeFromUnit2Start-timeFromUnit1Start;

                    timeDifferencesPerCycle(unit1Num,unit2Num,ci)=unit2TimingCurrCycle-unit1TimingCurrCycle;
                    
                    
                    %{
                    if(unit1Num~=unit2Num && abs(timeDifferencesPerCycle(unit1Num,unit2Num,ci))<1e-10)
                        disp('')
                        
                    end
                    %}
                    
                    if(~isnan(timeDifferencesPerCycle(unit1Num,unit2Num,ci)))
                        numCyclesCofiring(unit1Num,unit2Num,di)=numCyclesCofiring(unit1Num,unit2Num,di)+1;
                    end
                end
                
                alreadyProcessed(unit1Num,unit2Num)=true;
                alreadyProcessed(unit2Num,unit1Num)=true;
                
            end
            unit1Num
            %numCyclesCofiring(unit1Num,unit1Num,di)=NaN;
            numCyclesMaxPerRow(unit1Num)=numCyclesCofiring(unit1Num,unit1Num,di);
            %numCyclesThreshPerRow(unit1Num)=0.1*numCyclesCofiring(unit1Num,unit1Num,di);
            %numCyclesThreshPerRow(unit1Num)=0.3*numCyclesCofiring(unit1Num,unit1Num,di);
            %numCyclesThreshPerRow(unit1Num)=0.2*numCyclesCofiring(unit1Num,unit1Num,di);
            %numCyclesThreshPerRow(unit1Num)=0.15*numCyclesCofiring(unit1Num,unit1Num,di);
            numCyclesThreshPerRow(unit1Num)=cofiringThresh*numCyclesCofiring(unit1Num,unit1Num,di);
        end
        
        earlyThetaTimeDiffs=[];
       lateThetaTimeDiffs=[];
       numLatePhasePairs=0;
                    numEarlyPhasePairs=0;
                    
       allOverlappingPairCenterDiffs=[];
       
       speedPerCycle=currThetaSeqData.speedPerCycle;
       timeOfCycleTrough=currThetaSeqData.timeOfCycleTrough;
       
        for unit1Num=1:numUnits
            for unit2Num=1:numUnits %vs upper triangle matrix  
                
                %not overlapping
                if(numCyclesCofiring(unit1Num,unit2Num,di) < numCyclesThreshPerRow(unit1Num))
                    continue
                end
                
                center1=nanmedian(currFieldCentersPerCycle(:,unit1Num));
                center2=nanmedian(currFieldCentersPerCycle(:,unit2Num));
                 currSpaceDiffMag=abs(center2-center1);
                 
                 %allOverlappingPairCenterDiffs=[allOverlappingPairCenterDiffs;currSpaceDiffMag]; 
                 
                 if(currSpaceDiffMag>maxSpaceDiff*spaceFracMax || currSpaceDiffMag<maxSpaceDiff*spaceFracMin)
                     continue
                 end
                 
                 if(~isnan(currSpaceDiffMag))
                   allOverlappingPairCenterDiffs=[allOverlappingPairCenterDiffs;currSpaceDiffMag]; 
                 end
                 %{
                  if(currSpaceDiffMag<maxSpaceDiff*spaceFracMin)
                        continue
                  end
                 %}
                 
                   thisCellEarlyPhaseDiffs=[];
                    thisCellLatePhaseDiffs=[];
                    
                    thisCellPhaseDiffPerCycle=NaN(numCycles,1);
                    %{
                    figure(101)
                    close(figure(101))
                    %}
                    
                    
                    unit1SpikeTimePerCycle=NaN(numCycles,1);
                    unit2SpikeTimePerCycle=NaN(numCycles,1);
                    
                    %avgLeadingFieldTiming
                for ci=1:numCycles
                    
                    if(speedPerCycle(ci)<minSpeedForThetaDiff || speedPerCycle(ci)>maxSpeedForThetaDiff)
                        continue
                    end
                    
                    unit1TimingCurrCycle=currThetaSeqData.avgSpikeTimesByCycleAndUnit(ci,unit1Num,di);
                    unit2TimingCurrCycle=currThetaSeqData.avgSpikeTimesByCycleAndUnit(ci,unit2Num,di);
                    
                    currCycleStartTime=currThetaSeqData.timeOfCycleTrough(ci);
                    
                    %unit1SpikeTimePerCycle(ci)=unit1TimingCurrCycle+currCycleStartTime;
                    %unit2SpikeTimePerCycle(ci)=unit2SpikeTimePerCycle+currCycleStartTime;
                    
                    currPhaseDiffMag=abs(timeDifferencesPerCycle(unit1Num,unit2Num,ci));
                    
                    thisCellPhaseDiffPerCycle(ci)=currPhaseDiffMag;
                    
                    if(unit1Num~=unit2Num && unit1TimingCurrCycle > lateCycleTimeThresh && unit2TimingCurrCycle >lateCycleTimeThresh)   
                        %lateThetaTimeDiffs=[lateThetaTimeDiffs; currPhaseDiffMag/currSpaceDiffMag];
                        %lateThetaTimeDiffs=[lateThetaTimeDiffs; currPhaseDiffMag];
                        thisCellLatePhaseDiffs=[thisCellLatePhaseDiffs; currPhaseDiffMag];
                        
                        if(isnan(currPhaseDiffMag))
                            %continue
                        end
                        
                    end
                    
                    if(unit1Num~=unit2Num &&  unit1TimingCurrCycle < earlyCycleTimeThresh && unit2TimingCurrCycle <earlyCycleTimeThresh )
                        %earlyThetaTimeDiffs=[earlyThetaTimeDiffs; currPhaseDiffMag/currSpaceDiffMag];
                        %earlyThetaTimeDiffs=[earlyThetaTimeDiffs; currPhaseDiffMag];
                        thisCellEarlyPhaseDiffs=[thisCellEarlyPhaseDiffs; currPhaseDiffMag];
                        if(isnan(currPhaseDiffMag))
                            %continue
                        end
                       
                    end
                    %drawnow
                end 
                
                %{
                figure(112)
                plotRasterStyle(eventTimes,rowNum,startTick,endTick,tickColor)
                hold on
                plotRasterStyle(unit2SpikeTimePerCycle,rowNum,startTick,endTick,tickColor)
                %}
                
                %unit1EntryTimePerCycle=squeeze(currThetaSeqData.unitEntryTimePerCycleAndField(:,unit1Num,di));
                %{
                figure(111)
                plot(timeOfCycleTrough-unit1EntryTimePerCycle,thisCellPhaseDiffPerCycle)
                hold on
                xlim([0 3])
                drawnow
                %}
                earlyThetaTimeDiffs=[earlyThetaTimeDiffs; nanmean(thisCellEarlyPhaseDiffs)];
                lateThetaTimeDiffs=[lateThetaTimeDiffs; nanmean(thisCellLatePhaseDiffs)];
                
                %{
                length(earlyThetaTimeDiffs)
                length(lateThetaTimeDiffs)
                
                if(length(lateThetaTimeDiffs)>175)
                disp('')
                end
                %}
                
                
                if(isnan(nanmean(thisCellEarlyPhaseDiffs)) || isnan(nanmean(thisCellLatePhaseDiffs)))
                    continue
                end
                currRasterFig=figure(131+di)
                
                earlyDiff=nanmean(thisCellEarlyPhaseDiffs)*1000;
                lateDiff=nanmean(thisCellLatePhaseDiffs)*1000;
                
                if(exist('allPlottedLateDiffsSortedSaved','var'))
                    
                    %plotIdx=find(abs(lateDiff-allPlottedLateDiffsSortedSaved)<0.001);
                     %plotIdx=find(abs(abs(earlyDiff-lateDiff)-allPlottedLateDiffsSortedSaved)<0.001);
                     plotIdx=find(abs((earlyDiff+lateDiff)/2-allPlottedLateDiffsSortedSaved)<0.0001);

                     
                else
                    plotIdx=numLatePhasePairs;
                end

                %subplot(2,1,2)
                        %plotRasterStyle(0,length(lateThetaTimeDiffs),NaN,NaN,NaN,3)
                     
                         numLatePhasePairs=numLatePhasePairs+1;
                        if(~plotVert)
                         pLate=plotRasterStyle(lateDiff,numLatePhasePairs,NaN,NaN,[1 0 0],3);
                        else
                            pLate=plotRasterStyleVert(lateDiff,plotIdx,NaN,NaN,[1 0 0],3)
                        end 
                        hold on
                        xlim([0 cycleDur*lastFrac*0.6]*1000)
                       
                        ylim([0 numLatePhasePairs])
                        
                 %subplot(2,1,1)
                        %plotRasterStyle(0,length(earlyThetaTimeDiffs),NaN,NaN,NaN,1)
                        %hold on
                        %pEarly=plotRasterStyle(earlyDiff,numEarlyPhasePairs,NaN,NaN,NaN,3)
                            numEarlyPhasePairs=numEarlyPhasePairs+1;
                        if(plotVert)
                            pEarly=plotRasterStyleVert(earlyDiff,plotIdx,NaN,NaN,NaN,3)
                        else
                            pEarly=plotRasterStyle(earlyDiff,numEarlyPhasePairs,NaN,NaN,NaN,3)

                        end
                      
                    
                        
                         if(~plotVert)
                             xlim([0 cycleDur*lastFrac*0.6]*1000)
                            ylim([0 numEarlyPhasePairs])
                            xlabel('Avg spike time difference (msec)')
                            ylabel('Cell pair no.')
                         else
                            ylim([0 cycleDur*lastFrac*0.6]*1000)
                            xlim([0 numEarlyPhasePairs])
                            ylabel('Avg spike time difference (msec)')
                            xlabel('Cell pair no.')
                         end
                        
                        if(earlyDiff > lateDiff)
                            lineColor='b';
                        else
                            lineColor='g';
                        end
                        
                        if(~plotVert)
                            plot([earlyDiff lateDiff],[numEarlyPhasePairs numEarlyPhasePairs]-1,lineColor,'LineWidth',2)
                        else
                            plot([plotIdx plotIdx],[earlyDiff lateDiff],lineColor,'LineWidth',2)
                        end
                        
                        %allPlottedLateDiffs=[allPlottedLateDiffs; abs(earlyDiff-lateDiff)];
                        allPlottedLateDiffs=[allPlottedLateDiffs; (earlyDiff+lateDiff)/2];

                %{
                figure(101); hEarly=histogram(thisCellEarlyPhaseDiffs,100)
                hold on
                hLate=histogram(thisCellLatePhaseDiffs,100)
                hEarly.FaceColor='k';
                hLate.FaceColor='r';
                drawnow
                pause 0.1
                %}
            end %unit 2 loop
        end %unit 1 loop
        [allPlottedLateDiffsSorted,sortByLateDiffIdx]=sort(allPlottedLateDiffs);
        save(sprintf('plottedLateDiffs%s_%d',currSessName,di),'allPlottedLateDiffsSorted','sortByLateDiffIdx')
        
       setFigFontTo(18)
       try
       legend([pEarly pLate],'early phase','late phase')
       legend boxoff
       
       if(~plotVert)
      xlim([0 36.5])
       else
            ylim([0 36.5])
       end
      
      box off
        close(currRasterFig)
         end
        %figure; histogram(allOverlappingPairCenterDiffs,50)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %plot time differences vs speed, pos, time for each unique pair
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        speedColors=magma(100);
        speedColors=parula(100);
        
        %edges=linspace(0,0.04,30); 
        %edges=linspace(0,0.05,100);
         %edges=linspace(0,0.05,numBins);
         edges=linspace(0,cycleDur*lastFrac*0.6,numBins);
         %edges=linspace(0,cycleDur*lastFrac*0.55,numBins);
         
         nonNaNIdxes=~isnan(earlyThetaTimeDiffs) & ~isnan(lateThetaTimeDiffs);
         
         numEarlyPairs=sum(~isnan(earlyThetaTimeDiffs(nonNaNIdxes)));
         numLatePairs=sum(~isnan(lateThetaTimeDiffs(nonNaNIdxes)));
         
         if(~isempty(earlyThetaTimeDiffs(nonNaNIdxes)))
            [p,h,stats] = signrank(earlyThetaTimeDiffs(nonNaNIdxes),lateThetaTimeDiffs(nonNaNIdxes));
            p
            if(isfield(stats,'zval'))
            zval=stats.zval
            signedrank=stats.signedrank
            end
            %{
            Ach10_25_13 rightward
            p = 4.5637e-04
              zval: 3.5051
             signedrank: 529
            
            Ach10_25_13 leftward
            p =0.0074
          zval: 2.6795
    signedrank: 281
            
            Ach11 leftward
            p =5.5166e-15
            zval =7.8145
            signedrank = 4225
            
            %}
         end
         
         masterEarlyDiffs=[masterEarlyDiffs; earlyThetaTimeDiffs(:)];
         masterLateDiffs=[masterLateDiffs; lateThetaTimeDiffs(:)];
         
        figure;Ne=histcounts((earlyThetaTimeDiffs),edges); 
        plot(1000*edgesToBins(edges),smooth(Ne/sum(Ne),5),'k-','LineWidth',5);
        hold on; 
        Nl=histcounts((lateThetaTimeDiffs),edges); 
        plot(1000*edgesToBins(edges),smooth(Nl/sum(Nl),5),'r-','LineWidth',4);
        xlabel('Avg spike time difference in theta cycle (msec)')
        ylabel('Probability')
        legend(sprintf('first %d%% of cycle (n=%d pairs)',round(lastFrac*100),numEarlyPairs),sprintf('last %d%% of cycle (n=%d pairs)',round(lastFrac*100),numLatePairs))
        title(removeUnderscores(sprintf('%s, overlapping field pairs < %d cm apart across all theta cycles, %s',currSessName,round(spaceFracMax*100),currDiStr)))
        setFigFontTo(18)
        legend boxoff
        box off
        saveas(gcf,sprintf('%s_%sLateVsEarlyCyclePhaseDiffs.tif',currSessName,currDiStr))
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %just late vs early cycle
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if(si==length(sessionNames) && di==2)
            figure;Ne=histcounts((masterEarlyDiffs),edges); 
            plot(1000*edgesToBins(edges),smooth(Ne/sum(Ne),4),'k-','LineWidth',5);
            hold on; 
            Nl=histcounts((masterLateDiffs),edges); 
            plot(1000*edgesToBins(edges),smooth(Nl/sum(Nl),5),'r-','LineWidth',4);
            
            nonNaNIdxesMaster=~isnan(masterEarlyDiffs) & ~isnan(masterLateDiffs);
            
            [p,h,stats] = signrank(masterEarlyDiffs(nonNaNIdxesMaster),masterLateDiffs(nonNaNIdxesMaster))
            %{
            p =4.6025e-18
            zval: 8.6628
            signedrank: 14069
            nanmean(masterEarlyDiffs(nonNaNIdxesMaster))*1000 =

                   18.4912 msec

           nanmean(masterLateDiffs(nonNaNIdxesMaster))*1000 =

                   13.7497 msec
            %}
            xlabel('Avg spike time difference in theta cycle (msec)')
            ylabel('Probability')
            legend(sprintf('first %d%% of cycle (n=%d pairs)',round(lastFrac*100),length(masterEarlyDiffs(nonNaNIdxesMaster))),sprintf('last %d%% of cycle (n=%d pairs)',round(lastFrac*100),length(masterLateDiffs(nonNaNIdxesMaster))))
            title(removeUnderscores(sprintf('All overlapping field pairs < %d cm apart across all theta cycles',round(spaceFracMax*100))))
            setFigFontTo(18)
            legend boxoff
            saveas(gcf,sprintf('AllLateVsEarlyCyclePhaseDiffs.tif'))
        end
        
       
        %continue
        
        figure
       
        for ui=1:numUnits
            for uii=(ui+1):numUnits
                currPairTimeDiffsPerCycle=squeeze(timeDifferencesPerCycle(ui,uii,:));
                
                currPairTimeInMetaFieldPerCycle=squeeze(currTimeInMetaFieldPerCycle(ui,uii,:));
                
                center1PerCycle=currFieldCentersPerCycle(:,ui);
                center2PerCycle=currFieldCentersPerCycle(:,uii);
                
                if(numCyclesCofiring(ui,uii,di)<numCyclesThreshPerRow(ui) || numCyclesCofiring(ui,uii,di)<numCyclesThreshPerRow(uii))
                    continue
                end
                
                currPairSpaceDiffPerCycle=center2PerCycle-center1PerCycle;
                
                if(isnan(max(currPairTimeDiffsPerCycle)))
                    continue
                end
                
                validIdxes=~isnan(currPairTimeDiffsPerCycle) & abs(speedPerCycle)>minSpeed;
                currSpeedPerCycle=speedPerCycle(validIdxes);
                currPairTimeDiffsPerCycle=currPairTimeDiffsPerCycle(validIdxes);
                currPairSpaceDiffPerCycle=currPairSpaceDiffPerCycle(validIdxes);
                currPairTimeInMetaFieldPerCycle=currPairTimeInMetaFieldPerCycle(validIdxes);
                
                thisUnitMeanSpeed=nanmean((currSpeedPerCycle));
                speedColorIdx=ceil(thisUnitMeanSpeed*100);
                thisUnitMeanTimeDiff=nanmean(currPairTimeDiffsPerCycle);
                thisUnitMeanSpaceDiff=nanmean(currPairSpaceDiffPerCycle);
                
                %subplot(2,1,1)
                if(zScorePlot)
                    %subplot(2,2,1)
                      subplot(2,1,1)
                    plot(zscoreLFP(currSpeedPerCycle),zscoreLFP(currPairTimeDiffsPerCycle),'k.')
                    xlabel('Running speed (Z)')
                    ylabel('Time difference within theta cycle (Z)')
                    xlim([-3 3])
                     ylim([-3 3])
                     axis square
                      hold on
                      
                     nonNanIdxDiffs=~isnan(currSpeedPerCycle) & ~isnan(currPairTimeDiffsPerCycle);
                     rMat=corrcoef(currSpeedPerCycle(nonNanIdxDiffs),currPairTimeDiffsPerCycle(nonNanIdxDiffs));
                     speedVsTimeDiffR=rMat(1,2);
                     
                     allSpeedVsTimeDiffR=[allSpeedVsTimeDiffR; speedVsTimeDiffR];
                     %{
                     subplot(2,2,2)
                     histogram(allSpeedVsTimeDiffR,30)
                     
                     %}
                else
                    subplot(2,1,1)
                    plot(nanmean(currSpeedPerCycle),nanmean(currPairTimeDiffsPerCycle),'k.','MarkerSize',20)
                    xlabel('Avg running speed (m/s)')
                    ylabel('Avg time difference within theta cycle (sec)')
                     xlim([0.2 1])
                     if(contains(currSessName,'11012013'))
                        ylim([-0.04 0.03])
                     else
                          ylim([-0.06 0.06])
                     end
                      hold on
                 %ylim([-Inf Inf])     
                end
               
                
                subplot(2,1,2)
                 %plot(currPairSpaceDiffPerCycle,currPairTimeDiffsPerCycle,'k.')
                 %{
                 plot(thisUnitMeanSpaceDiff,thisUnitMeanTimeDiff,'.','MarkerSize',20,'Color',speedColors(speedColorIdx,:))
                 colorbar
                 caxis([0 1])
                 %}
                %scatter(thisUnitMeanSpaceDiff,thisUnitMeanTimeDiff,20,thisUnitMeanSpeed
                 
                 plot(thisUnitMeanSpaceDiff,thisUnitMeanTimeDiff,'k.','MarkerSize',20)
                    xlabel('Spatial difference between field centers (m)')
                    ylabel('Avg time difference within theta cycle (sec)')
                    
                hold on
                xlim([-0.5 0.5])
                     if(contains(currSessName,'11012013'))
                        ylim([-0.05 0.05])
                     else
                          ylim([-0.06 0.06])
                     end
                 %ylim([-Inf Inf])
                 
                 %{
                 subplot(3,1,3)
                  plot(nanmean(currPairTimeInMetaFieldPerCycle),nanmean(currPairTimeDiffsPerCycle),'k.','MarkerSize',20)
                  

                  xlabel('Time since field entry (sec)')
                    ylabel('Time difference within theta cycle (sec)')
                   xlim([-1 1])
                       %ylim([-0.06 0.06])
                     %axis square
                     hold on
                 %}
                 
                 
            end
        end
        uberTitle(removeUnderscores(sprintf('%s, all overlapping field pairs across all theta cycles, %s',currSessName,currDiStr)))
                %maxFigFourthWidthFullHeight
                maxFigHalfWidth
                setFigFontTo(18)
                saveas(gcf,sprintf('%s_%sSpeedInvariantSpatialSeqCoding.tif',currSessName,currDiStr))
                
               
                %speedDiffRedges=linspace(-0.5,0.5,50);
                speedDiffRedges=linspace(-1,1,100);
                speedDiffBinCenters=edgesToBins(speedDiffRedges);
                %pSpeedThetaDiff=getProbDist(allSpeedVsTimeDiffR,speedDiffRedges,1,0);
                 pSpeedThetaDiff=getProbDist(allSpeedVsTimeDiffR,speedDiffRedges,0,0);
                
                
                 figure
                plot(speedDiffBinCenters,pSpeedThetaDiff,'k-')
                xlabel('Speed-phase difference correlation coefficient')
                legend(sprintf('n=%d field pairs',length(allSpeedVsTimeDiffR)))
                axis tight
        
                hold on; plot([0 0],ylim,'k--')
                ylabel('Probability')
                setFigFontTo(18)
                xlim([-1 1])
                        disp('')
                        box off
                        legend boxoff
                        maxFigHalfHalfWidth
                
                
                %{
                nanmean(allSpeedVsTimeDiffR) = -0.0083
                getSEMacrossRows(allSpeedVsTimeDiffR) = 0.0085
                
                Signed rank test different from 0:
                p =0.3567
                zval: -0.9217
                signedrank: 21189
                
                %}
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %plot pairwise time differences heat map
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %{
        avgTimeDiffPerCycle=nanmean(timeDifferencesPerCycle,3);
        
        numCyclesCofiringThisDir=numCyclesCofiring(:,:,di);
        
        for u=1:numUnits
            noOverlapThisRow=numCyclesCofiringThisDir(u,:)<numCyclesThreshPerRow(u);
            numCyclesCofiringThisDir(u,:)=numCyclesCofiringThisDir(u,:)/ numCyclesMaxPerRow(u);
            numCyclesCofiringThisDir(u,noOverlapThisRow)=NaN;
            %numCyclesCofiringThisDir(noOverlapThisRow,u)=NaN;
            
            avgTimeDiffPerCycle(u,noOverlapThisRow)=NaN;
            %avgTimeDiffPerCycle(noOverlapThisRow,u)=NaN;
            
        end
        
        avgTimeDiffPerCycle=avgTimeDiffPerCycle(sortedMeanCenterIdxPerDi(:,di),sortedMeanCenterIdxPerDi(:,di));
        numCyclesCofiringThisDir=numCyclesCofiringThisDir(sortedMeanCenterIdxPerDi(:,di),sortedMeanCenterIdxPerDi(:,di));
        
        removeRowIdxes=[];
        removeColIdxes=[];
        for n=1:numUnits
            if(isnan(max(avgTimeDiffPerCycle(n,:))) || sum(~isnan(avgTimeDiffPerCycle(n,:)))==1)
                removeRowIdxes=[removeRowIdxes; n];
            end
            
            if(isnan(max(avgTimeDiffPerCycle(:,n))) || sum(~isnan(avgTimeDiffPerCycle(:,n)))==1)
                removeColIdxes=[removeColIdxes; n];
            end
        end
        
        
        avgTimeDiffPerCycle(removeRowIdxes,:)=[];
        avgTimeDiffPerCycle(:,removeColIdxes)=[];
        numCyclesCofiringThisDir(removeRowIdxes,:)=[];
        numCyclesCofiringThisDir(:,removeColIdxes)=[];
        
        figure; 
        subplot(1,2,1)
        imagesc(avgTimeDiffPerCycle)
        colormap(jet)
        colorbar
        caxis([-0.065 0.065])
        
        subplot(1,2,2)
        colormap(gca,parula)
        imagesc(numCyclesCofiringThisDir)
        colorbar
        %}
        
    end   
end
        
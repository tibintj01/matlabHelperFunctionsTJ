close all;  clc
%dataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles';

clearvars -except spikeDataPerField

offsetAndDiffVsSpeedDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3Images/offsetAndDiffVsSpeedIndThetaSeqMin2FieldsCntrlDist';

exampleThetaSeqVsSpeedDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3Images/exampleThetaSeqVsSpeedNoRangeRestr';
touchDir(exampleThetaSeqVsSpeedDir)

%nonCircStats=1;
nonCircStats=0;

usePeakPhaseAsSummary=0;

controlForDistanceAcrossLaps=1;

touchDir(offsetAndDiffVsSpeedDir)
tic
disp('loading spike data per field')
if(~exist('spikeDataPerField','var'))
    spikeDataPerField=load('allUnitInFieldGoodSpeedSpikePhasesAndTimes.mat');
end
toc

maxRunSpeed=1;

minSpeed=0.1;
minSpeed=0.125;
%maxSpeed=0.7;
maxSpeed=0.7;

speedInFieldsPerSec=1;

if(speedInFieldsPerSec)
    %minSpeed=0.2;
    %minSpeed=0.3;
    minSpeed=0.25;
    %maxSpeed=1.2;
     %maxSpeed=1.25;
     maxSpeed=1.2;
end

%zScoreSpeed=1;
zScoreSpeed=0;

   minZspeed=-2.5;
    maxZspeed=2;
    
processedDataDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData';
unitDataDir=fullfile(processedDataDir,'unitSpikeInfoStructs');
dataDir=unitDataDir;

saveThetaSeqImgDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3Images/individualThetaSequencesSimilarSizeMinTrajLengthHalfWidth';
touchDir(saveThetaSeqImgDir)

posInfoDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData/positionTimeInfos';
perCycleDataDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData/perFieldPerCycleData_Dec16_2021';

thetaSequenceImgDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData/thetaSequenceImgDir';
touchDir(thetaSequenceImgDir)


%cycleDataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/';


setTightSubplots
sqlDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/sqlQueries';

sqlSessionFiles=getRegexFileNames(sqlDir,'ec*csv');

numSessions=length(sqlSessionFiles);
sessionNames=cell(numSessions,1);

for si=1:numSessions
    sessionNames{si}=getSubstrBtwnSubstrs(sqlSessionFiles{si},'','_sql');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%loop through each session and load per cycle spike time data across units
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

startSi=1;
useCircMeasureSpread=1;
useCircMeasureSpread=0;

if(zScoreSpeed)
    minSpeedStdDevAcrossLaps=0.05;
else
    minSpeedStdDevAcrossLaps=0;
end
%maxSpeedStdDevWithinTrav=0.1;
maxSpeedStdDevWithinTrav=1;

minNumPtsInSpeedCategory=2;
minNumPtsInSpeedCategory=1;

%minSpeedStdDevAcrossLaps=0.15;
minNumOccurences=5;
minNumFieldsInSeq=2;
%minNumFieldsInSeq=4;
showPlots=1;
%showPlots=0;

setTightSubplots_Spacious

fH=figure;
masterTravSpeed=[];
masterPhaseDiff=[];
masterMeanPhase=[];
masterPhaseDiffSpeedRhos=[];
masterMeanPhaseSpeedRho=[];
masterPhaseDiffSpeedSlope=[];
masterMeanPhaseSpeedSlope=[];

masterPhaseDiffHighMinusLowSpeed=[];
masterMeanPhaseHighMinusLowSpeed=[];
masterMeanDistFracHighMinusLowSpeed=[];
indFieldHighSpeedPrctileCutoff=80;
indFieldLowSpeedPrctileCutoff=20;

maxAllowableSlopeMagnitude=720;
maxAllowableSlopeMagnitude=Inf;


startSi=1;

for si=startSi:length(sessionNames) 
    si
    currSessName=sessionNames{si};
    currSeqDataFilePath=fullfile(perCycleDataDir,sprintf('perFieldPerCycleTimingData%s.mat',currSessName));
    
    %missing seq files?
    if(exist(currSeqDataFilePath,'file'))
        currThetaSeqData=load(currSeqDataFilePath);
    else
        continue
    end
    
    %{
    if(~strcmp(currSessName,'ec013.626'))
        continue
    end
    %}
    
    
    if(isempty(currThetaSeqData.thetaSequenceInfoPerSequence))
        continue
    end
    
    %{
    if(~isfield(currThetaSeqData,'thetaSequenceInfoPerSequenceAtLeast4Fields'))
        continue
    end
    %}
    
    thetaSequenceInfoPerSequence=currThetaSeqData.thetaSequenceInfoPerSequence;
       % thetaSequenceInfoPerSequence=currThetaSeqData.thetaSequenceInfoPerSequenceAtLeast4Fields;

    
    seqStrs=fieldnames(thetaSequenceInfoPerSequence);
    
    for ti=1:length(seqStrs)
        currSeqData=thetaSequenceInfoPerSequence.(seqStrs{ti});
        numOccurences=currSeqData.numOccurences;
       
        avgFieldWidth=currSeqData.avgFieldWidth;
        numFieldsInSeq=length(currSeqData.fieldCenterSeq);
        travSpeedPerOcc=cell2num(currSeqData.travAvgSpeedPerOcc);
        
        seqDirStr=currSeqData.seqDirStr;
        
        [sortedSpeedPerOcc,sortBySpeedOccIdxes]=sort(travSpeedPerOcc);
        
      
        allPhaseSeqPerSpeedSortedOcc=currSeqData.allPhaseSeqPerOcc(sortBySpeedOccIdxes);
        meanPhaseSeqPerSpeedSortedOcc=currSeqData.meanPhaseSeqPerOcc(sortBySpeedOccIdxes);
        
        allPhaseSeqPerOriginalOcc=currSeqData.allPhaseSeqPerOcc;
        
        fieldMaxRatePosPerField=NaN(numFieldsInSeq,1);

        for ai=1:numFieldsInSeq
            currFieldFiringRatePerPos=currSeqData.firingRatePerPosPerFieldInSeq{ai};

            [~,currFieldMaxRatePosBin]=max(smooth(currFieldFiringRatePerPos(:,2)));

            currFieldMaxRatePos=currFieldFiringRatePerPos(currFieldMaxRatePosBin,1);

            fieldMaxRatePosPerField(ai)=currFieldMaxRatePos;
        end
        
        if(strcmp(seqDirStr,'rightward'))
            [~,fieldIdxToCenterRank]=sort(fieldMaxRatePosPerField);
        else(strcmp(seqDirStr,'leftward'))
            [~,fieldIdxToCenterRank]=sort(fieldMaxRatePosPerField,'descend');
        end
                
        
        
        if(numFieldsInSeq<minNumFieldsInSeq)
            continue
        end
        
        
        
        speedStdDevAcrossLaps=nanstd(cell2num(currSeqData.travAvgSpeedPerOcc));
        
        if(numOccurences>=minNumOccurences && speedStdDevAcrossLaps>=minSpeedStdDevAcrossLaps...
                & numFieldsInSeq>=minNumFieldsInSeq)
            %currSeqData

            %fH=figure; 
        distFracColors=jet(numFieldsInSeq);
        allFieldDistFracs=[currSeqData.fieldDistFracSeqPerOcc{:}];
        
        repFracTol=0.05;
        %repFracTol=0.2/3;
        %repFracTol=0.075;
        %repFracTol=0.025;
        %repFracTol=0.1;
        
        isControlledDistPerOcc=ones(numOccurences,1);
        for fi=1:numFieldsInSeq
            currFieldDistFracs=allFieldDistFracs(fi,:);
            
            gaussKernelSD=0.2/3; %3 of these is approximately full width
            %gaussKernelSD=0.15/3; %3 of these is approximately full width
            %gaussKernelSD=0.3/3; %3 of these is approximately full width
            distrEdges=linspace(0,1,100);
            [currPdist,currFieldRepFrac] = getKernelDistrAndPeakVal(currFieldDistFracs,gaussKernelSD,distrEdges);
            
            currFieldCloseToRepFracOccIdxes=currFieldDistFracs>=currFieldRepFrac-repFracTol & currFieldDistFracs<=currFieldRepFrac+repFracTol;
            
            isControlledDistPerOcc=isControlledDistPerOcc(:) & (currFieldCloseToRepFracOccIdxes(:)); %must be good for all fields
            
            currFieldDistFracs(~isControlledDistPerOcc)=NaN;
            
            %{
            hold on
            plot(1:length(currFieldDistFracs),currFieldDistFracs,'.','Color',distFracColors(fi,:),'MarkerSize',30)
             hold on
             hline(currFieldRepFrac,'k--',3)
            %}
            
            
        end
        meanDistFracPerField=nanmean([currSeqData.fieldDistFracSeqPerOcc{:}],2);
       
        %isControlledDistPerOcc
        %hline(meanDistFracPerField,'k-',3)
        
        %close all
        %continue
        
        
        
        if(usePeakPhaseAsSummary)
            peakPhaseSeqPerOcc=currSeqData.peakPhaseSeqPerOcc;
        else
            
              peakPhaseSeqPerOcc=currSeqData.meanPhaseSeqPerOcc;
        end

            %travAvgSpeedPerOcc=currSeqData.travAvgSpeedPerOcc;
            %travSpeedPerOcc=cell2num(travAvgSpeedPerOcc);
            
            if(speedInFieldsPerSec)
                travSpeedPerOcc=travSpeedPerOcc/avgFieldWidth; %from m/s to field widths per second
            end
            
            if(speedInFieldsPerSec)
                sortedSpeedPerOcc=sortedSpeedPerOcc/avgFieldWidth;
            end
            travSpeedTracePerOcc=currSeqData.travSpeedTracePerOcc;
            numFieldsInSeq=length(currSeqData.fieldCenterSeq);
            
            avgPhaseDiffPerOcc=NaN(numOccurences,1);
            avgPhasePerOcc=NaN(numOccurences,1);
            speedStdDevWithinTrav=NaN(numOccurences,1);
            
            for oi=1:numOccurences
                
                if(~isControlledDistPerOcc(oi))
                    %continue
                end
                %for fi=1:numFieldsInSeq
                    currOccPhaseSeq=peakPhaseSeqPerOcc{oi};
                    
                    %avgPhaseDiffPerOcc(oi)=nanmean(angdiffDeg(currOccPhaseSeq)); %measure of phase spread
                   if(useCircMeasureSpread)
                         avgPhaseDiffPerOcc(oi)=1-circ_r(ang2rad(currOccPhaseSeq(:)));
                   else
                         avgPhaseDiffPerOcc(oi)=circMeanDegNoMod360(angdiffDeg(currOccPhaseSeq)); %measure of phase spread

                   end

                    if(nonCircStats)
                        avgPhasePerOcc(oi)=nanmean(currOccPhaseSeq);             %measure of central phase
                    else
                        avgPhasePerOcc(oi)=circMeanDeg(currOccPhaseSeq);             %measure of central phase
                    end

                     %avgPhasePerOcc(oi)=circMedianDeg(currOccPhaseSeq);             %measure of central phase
                     
                     speedStdDevWithinTrav(oi)=nanstd(travSpeedTracePerOcc{oi});

                %end
            end
            
            
            
            computeCorrIdxes=travSpeedPerOcc>=minSpeed & travSpeedPerOcc<=maxSpeed;
            
            constSpeedTravIdxes=speedStdDevWithinTrav <= maxSpeedStdDevWithinTrav;
            
            %Dec22
            %{
            avgPhaseDiffPerOcc(~constSpeedTravIdxes)=NaN;
            avgPhasePerOcc(~constSpeedTravIdxes)=NaN;
            %}
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %relative to ind cell speed cutoff
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            highSpeedMinCutoff=prctile(travSpeedPerOcc,indFieldHighSpeedPrctileCutoff);
            lowSpeedMaxCutoff=prctile(travSpeedPerOcc,indFieldLowSpeedPrctileCutoff);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %absolute cutoff
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            highSpeedMinCutoff=0.5;
            lowSpeedMaxCutoff=0.3;
            
             highSpeedMinCutoff=0.45;
            lowSpeedMaxCutoff=0.35;
            
            highSpeedMinCutoff=0.55;
            lowSpeedMaxCutoff=0.25;
            
             % highSpeedMinCutoff=0.4;
            %lowSpeedMaxCutoff=0.4;
            
             %highSpeedMinCutoff=0.5;
            %lowSpeedMaxCutoff=0.3;
            
             highSpeedMinCutoff=0.45;
            lowSpeedMaxCutoff=0.35;
             %highSpeedMinCutoff=0.6;
            %lowSpeedMaxCutoff=0.2;
             highSpeedMinCutoff=0.5;
            lowSpeedMaxCutoff=0.3;
            
            highSpeedMinCutoff=0.55;
            lowSpeedMaxCutoff=0.25;
            
            if(speedInFieldsPerSec) 
                highSpeedMinCutoff=1.1;
                %lowSpeedMaxCutoff=0.5;
                %highSpeedMinCutoff=0.8;
                %lowSpeedMaxCutoff=0.7;
                %lowSpeedMaxCutoff=0.8;
                lowSpeedMaxCutoff=0.75;
                
                highSpeedMinCutoff=1;
                lowSpeedMaxCutoff=0.4;
                
                 highSpeedMinCutoff=1;
                lowSpeedMaxCutoff=0.5;
                
                highSpeedMinCutoff=1.1; %p ***
                lowSpeedMaxCutoff=0.55;
                
                 highSpeedMinCutoff=1.05; %p **
                lowSpeedMaxCutoff=0.55;
                
                %highSpeedMinCutoff=1.05; %p ***
                %lowSpeedMaxCutoff=0.6;
            end
            
            %highSpeedMinCutoff=0.6;
            %lowSpeedMaxCutoff=0.2;
            %highSpeedIdxes=computeCorrIdxes & travSpeedPerOcc>=highSpeedMinCutoff;
            %lowSpeedIdxes=computeCorrIdxes & travSpeedPerOcc<=lowSpeedMaxCutoff;
            
             %highSpeedIdxes=travSpeedPerOcc>=highSpeedMinCutoff;
            %lowSpeedIdxes=travSpeedPerOcc<=lowSpeedMaxCutoff;
            
             highSpeedIdxes=constSpeedTravIdxes(:) & travSpeedPerOcc(:)>=highSpeedMinCutoff;
            lowSpeedIdxes=constSpeedTravIdxes(:) & travSpeedPerOcc(:)<=lowSpeedMaxCutoff;
            
            highSpeedDistFracs=currFieldDistFracs(highSpeedIdxes);
            lowSpeedDistFracs=currFieldDistFracs(lowSpeedIdxes);
            
            %figure place holder
            circFigH=[];
            
            enoughPtsInEachCategory=sum(highSpeedIdxes)>=minNumPtsInSpeedCategory && sum(lowSpeedIdxes)>=minNumPtsInSpeedCategory;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %plot summary of this field sequence vs speed
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %if(showPlots && enoughPtsInEachCategory)
            if(showPlots)
                showSpeedSortedRaster=1;
                figure
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 %rate fields
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                subplot(2,2,1)
                distColors=jet(numFieldsInSeq);
                
                
                
                
                for ai=1:numFieldsInSeq
                    currFieldFiringRatePerPos=currSeqData.firingRatePerPosPerFieldInSeq{ai};
                    
                   
                    plot(currFieldFiringRatePerPos(:,1),currFieldFiringRatePerPos(:,2),'Color',distColors(fieldIdxToCenterRank(ai),:),'LineWidth',3)
                    hold on
                end
                
                box off
                
                xlabel('Position (m)')
              ylabel('Firing rate (Hz)')
              axis tight
                
                colormap(jet(numFieldsInSeq))
               cb=colorbar;
               cb.Ticks=[];
                ylabel(cb,'Behavioral field order')
                caxis([1 numFieldsInSeq])
                
                title(sprintf('Place fields in behavioral sequence'))%, %s %s',currSessName,removeUnderscores(seqStrs{ti})))
                
               
                
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 %speed variability
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                subplot(2,2,3)
 
                speedEdges=linspace(minSpeed,maxSpeed,21);
                histogram(travSpeedPerOcc,speedEdges)
                hold on

                 pLh=vline(lowSpeedMaxCutoff,'k--',3);
                    pHh=vline(highSpeedMinCutoff,'k-',3);
                    
                    legend([pLh pHh],{sprintf('%.2f fields/sec',lowSpeedMaxCutoff),sprintf('%.2f fields/sec',highSpeedMinCutoff)},'Location','westoutside')
                    legend boxoff

                 %pS=vline(currMeanTravSpeed,'k--',2);
             
             %legend([pS],{'curr. speed'},'Location','northeast')%'eastoutside')
             %legend boxoff
            
             xlim([minSpeed,maxSpeed])
             
             if(speedInFieldsPerSec)
                 xlabel('traversal speed (fields/sec)')
             else
                xlabel('traversal speed (m/s)')
             end
            
            ylabel('traversal count')
             axis tight
             %axis square
             box off

             title('Speed distribution across laps')
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %raster of speed sorted theta sequences
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if(showSpeedSortedRaster)
                    
                    %subplot(1,2,2)
                    subplot(2,2,[2 4])

                    for oi=1:numOccurences

                        currOccAllPhaseSeq=allPhaseSeqPerSpeedSortedOcc{oi};
                        currOccMeanPhaseSeq=meanPhaseSeqPerSpeedSortedOcc{oi};

                        %rasterAlpha=0.3;
                        rasterAlpha=0.4;
                        for ofi=1:length(currOccAllPhaseSeq)
                            currOccFieldPhases=currOccAllPhaseSeq{ofi};
                            pRaster=plotRasterStyle(currOccFieldPhases,oi,NaN,NaN,[distColors(fieldIdxToCenterRank(ofi),:) rasterAlpha])
                            hold on
                            %pRaster.Color(4)=rasterAlpha;
                            %plot(currOccMeanPhaseSeq(ofi),oi,'.','Color',distColors(fieldIdxToCenterRank(ofi),:),'MarkerSize',round(20*20/numOccurences))
                            %plot(currOccMeanPhaseSeq(ofi),oi,'.','MarkerFaceColor',distColors(fieldIdxToCenterRank(ofi),:),'MarkerEdgeColor','k','MarkerSize',round(20*20/numOccurences))
                            %scatter(currOccMeanPhaseSeq(ofi),oi,round(50*20/numOccurences),'MarkerFaceColor',distColors(fieldIdxToCenterRank(ofi),:),'MarkerEdgeColor','k','LineWidth',1)
                            scatter(currOccMeanPhaseSeq(ofi),oi,round(200*20/numOccurences),'MarkerFaceColor',distColors(fieldIdxToCenterRank(ofi),:),'MarkerEdgeColor','k','LineWidth',1)

                            %scatter(x,y,sz,'MarkerEdgeColor',[0 .5 .5],...
                            %  'MarkerFaceColor',[0 .7 .7],...
                            %  'LineWidth',1.5)
                            %plot(currOccMeanPhaseSeq(ofi),oi,'ko','MarkerSize',round(100/numOccurences))

                        end
                    end
                    
                    if(lowSpeedMaxCutoff>=min(sortedSpeedPerOcc))
                        [~,lowSpeedMaxCutoffIdx]=min(abs(sortedSpeedPerOcc-lowSpeedMaxCutoff));
                        lowSpeedMaxCutoffIdx=round(lowSpeedMaxCutoffIdx)+0.5;
                    else
                        lowSpeedMaxCutoffIdx=0;
                    end
                    
                    if(highSpeedMinCutoff<=max(sortedSpeedPerOcc))
                        [~,highSpeedMinCutoffIdx]=min(abs(sortedSpeedPerOcc-highSpeedMinCutoff));
                        highSpeedMinCutoffIdx=round(highSpeedMinCutoffIdx)-0.5;
                    else
                        highSpeedMinCutoffIdx=length(sortedSpeedPerOcc)+0.5;
                    end
                    
                    pL=hline(lowSpeedMaxCutoffIdx,'k--',3);
                    pH=hline(highSpeedMinCutoffIdx,'k-',3);
                    
                    %legend([pL pH],{sprintf('%.2f fields/sec',lowSpeedMaxCutoff),sprintf('%.2f fields/sec',highSpeedMinCutoff)});%,'Location','eastoutside')
                    %legend boxoff
                    box off
                    xlim([0 360])
                    ylim([0 numOccurences+1])
                    
                      xlabel('Theta phase (deg)')
                    ylabel('Speed-sorted traversal no.')
                    
                        colormap(jet(numFieldsInSeq))
                   cb=colorbar;
                   cb.Ticks=[]
                    ylabel(cb,'Behavioral field order')
                    caxis([1 numFieldsInSeq])
                    %}

                else %show low speed vs high speed rasters



                    %allPhaseSeqPerOriginalOcc

                    highSpeedOccCount=0;
                    lowSpeedOccCount=0;

                    for oi=1:numOccurences

                        if(highSpeedIdxes(oi))
                            subplot(2,2,2)
                            currOccAllPhaseSeq=allPhaseSeqPerOriginalOcc{oi};
                            highSpeedOccCount=highSpeedOccCount+1;

                            for ofi=1:length(currOccAllPhaseSeq)
                                currOccFieldPhases=currOccAllPhaseSeq{ofi};
                                plotRasterStyle(currOccFieldPhases,highSpeedOccCount,NaN,NaN,distColors(fieldIdxToCenterRank(ofi),:))
                                hold on
                            end

                            xlim([0 360])
                            ylim([0 highSpeedOccCount+1])
                        end

                        if(lowSpeedIdxes(oi))
                            subplot(2,2,4)
                            currOccAllPhaseSeq=allPhaseSeqPerOriginalOcc{oi};
                            lowSpeedOccCount=lowSpeedOccCount+1;

                            for ofi=1:length(currOccAllPhaseSeq)
                                currOccFieldPhases=currOccAllPhaseSeq{ofi};
                                plotRasterStyle(currOccFieldPhases,lowSpeedOccCount,NaN,NaN,distColors(fieldIdxToCenterRank(ofi),:))
                                hold on
                            end

                            xlim([0 360])
                            ylim([0 lowSpeedOccCount+1])
                        end
                    end
                end
                title('All theta sequences vs running speed')
              

                
                
                 uberTitle(removeUnderscores(sprintf('%s, %s',currSessName,seqStrs{ti})))
                setFigFontTo(16)
                
                try
                    maxFig
                    saveas(gcf,fullfile(exampleThetaSeqVsSpeedDir,sprintf('thetaSequencesVsSpeed_%s_%s.png',currSessName,seqStrs{ti})))
                catch
                    maxFigMukkaalWidth
                    saveas(gcf,fullfile(exampleThetaSeqVsSpeedDir,sprintf('thetaSequencesVsSpeed_%s_%s.png',currSessName,seqStrs{ti})))

                end
                close all
                continue

            end

        
            
            if(showPlots && enoughPtsInEachCategory)
                subplot(1,2,1)
                %travSpeedPerOcc=cell2num(travAvgSpeedPerOcc);
                plot(travSpeedPerOcc,avgPhaseDiffPerOcc,'ko')

                hold on
                plot(travSpeedPerOcc(highSpeedIdxes),avgPhaseDiffPerOcc(highSpeedIdxes),'r.','MarkerSize',10)
                plot(travSpeedPerOcc(lowSpeedIdxes),avgPhaseDiffPerOcc(lowSpeedIdxes),'b.','MarkerSize',10)

                title('Avg phase diff in theta sequence vs traversal speed')

                 subplot(1,2,2)
                
                plot(travSpeedPerOcc,avgPhasePerOcc,'ko')
                hold on
                hold on
                plot(travSpeedPerOcc(highSpeedIdxes),avgPhasePerOcc(highSpeedIdxes),'r.','MarkerSize',10)
                plot(travSpeedPerOcc(lowSpeedIdxes),avgPhasePerOcc(lowSpeedIdxes),'b.','MarkerSize',10)
                
                title('Avg phase in theta sequence vs traversal speed')
                drawnow
                circFigH=fH;
            end
            
            if(false && showPlots && enoughPtsInEachCategory)
                subplot(1,2,1)
                box off
                 xlabel('Avg speed in traversal (m/s)')
                ylabel('Avg phase diff. in theta seq. (deg)')
                axis square
                
                %figure; 
                %close all
                 numImages=length(currSeqData.currThetaSeqSummaryImgPaths);
                 highSpeedImgCount=1;
                 lowSpeedImgCount=1;
                for ai=1:numImages
                    
                    if(highSpeedIdxes(ai))
                            figure(11)
                        subplot(floor(sqrt(sum(highSpeedIdxes))),floor(sqrt(sum(highSpeedIdxes)))+1,highSpeedImgCount)

                       imshow(currSeqData.currThetaSeqSummaryImgPaths{ai})
                       
                       highSpeedImgCount=highSpeedImgCount+1;
                   title(sprintf('Cycle %d',currSeqData.cycleNumPerOcc{ai}))
                    end
                    
                    if(lowSpeedIdxes(ai))
                            figure(12)
                        subplot(floor(sqrt(sum(lowSpeedIdxes))),floor(sqrt(sum(lowSpeedIdxes)))+1,lowSpeedImgCount)

                       imshow(currSeqData.currThetaSeqSummaryImgPaths{ai})
                      
                       
                       lowSpeedImgCount=lowSpeedImgCount+1;
                       title(sprintf('Cycle %d',currSeqData.cycleNumPerOcc{ai}))
                   
                    end
                    
                     
                    %autoArrangeFigures()
                    
                    
                end
                
                close all
            end
            
            
            
            circFigH=[]; %no circ corr lines
           
           %{ 
            [currSeqPhaseDiffSpeedRho,currSeqPhaseDiffSpeedPval,currSeqPhaseDiffSpeedSlopeDegPerMsec, ~]...
                =getCircCorrCoeff(travSpeedPerOcc(computeCorrIdxes),avgPhaseDiffPerOcc(computeCorrIdxes),circFigH);
            
            %}
            
            if(showPlots && enoughPtsInEachCategory)
                
                if(useCircMeasureSpread)
                    ylim([0 1])
                else
                    ylim([-180 180])
                end
            end
            
            
         
            if(showPlots && enoughPtsInEachCategory)
                subplot(1,2,2)
                box off
                xlabel('Avg speed in traversal (m/s)')
                ylabel('Avg phase of theta seq. (deg)')
                axis square
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %plot raster examples similar to Fig 1B
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if(numFieldsInSeq>=4)
                    %{
                    figure;
                    speedRank
                    
                    for oi=1:numOccurences
                        peakPhaseSeqPerOcc{oi}
                    end
                    %}
                    
                end
            end
            %{
            [currSeqMeanPhaseSpeedRho,currSeqMeanPhaseSpeedPval,currSeqMeanPhaseSpeedSlopeDegPerMsec, ~]...
                =getCircCorrCoeff(travSpeedPerOcc(computeCorrIdxes),avgPhasePerOcc(computeCorrIdxes),circFigH);
            %}
            
            
            
            
              if(showPlots && enoughPtsInEachCategory)
                ylim([0 360])
              end
            
              if(showPlots && enoughPtsInEachCategory)
                  uberTitle(removeUnderscores(sprintf('%s, %s',currSessName,seqStrs{ti})))
          
                  maxFigHalfHalfWidth
                  setFigFontTo(16)
                  saveas(gcf,fullfile(offsetAndDiffVsSpeedDir,sprintf('indivThetaSeqOffsetAndDiffVsSpeed_%s_Seq%s.png',currSessName,seqStrs{ti})))
                   cla
                  %close all
               
                  subplot(1,2,1)
                hold off
                subplot(1,2,2)
                hold off
               
              end
            
           
                
              %{
             if(abs(currSeqPhaseDiffSpeedSlopeDegPerMsec)<= maxAllowableSlopeMagnitude ...
                     && abs(currSeqMeanPhaseSpeedSlopeDegPerMsec)<=maxAllowableSlopeMagnitude)
                 
              %}
              %{
                 masterPhaseDiffSpeedRhos=[masterPhaseDiffSpeedRhos; currSeqPhaseDiffSpeedRho(:)];
                 masterMeanPhaseSpeedRho=[masterMeanPhaseSpeedRho; currSeqMeanPhaseSpeedRho(:)];
          
                 masterPhaseDiffSpeedSlope=[masterPhaseDiffSpeedSlope; currSeqPhaseDiffSpeedSlopeDegPerMsec(:)];
                 masterMeanPhaseSpeedSlope=[masterMeanPhaseSpeedSlope; currSeqMeanPhaseSpeedSlopeDegPerMsec(:)];
             %}

                 if(zScoreSpeed)
                   masterTravSpeed=[masterTravSpeed; zscoreLFP(travSpeedPerOcc(:))];
                 else
                    
                     masterTravSpeed=[masterTravSpeed; travSpeedPerOcc(:)];
                 end
                
                masterPhaseDiff=[masterPhaseDiff; avgPhaseDiffPerOcc(:)];
                masterMeanPhase=[masterMeanPhase; avgPhasePerOcc(:)];
            
             %end
             
             if(sum(highSpeedIdxes)<minNumPtsInSpeedCategory || sum(lowSpeedIdxes)<minNumPtsInSpeedCategory)
                continue
             end
            
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %get high speed and low speed values
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            meanDistFracHighSpeed=nanmean(highSpeedDistFracs);
             meanDistFracLowSpeed=nanmean(lowSpeedDistFracs);
            
            
            meanPhaseDiffHighSpeed = circMeanDeg(avgPhaseDiffPerOcc(highSpeedIdxes));
             meanPhaseDiffLowSpeed = circMeanDeg(avgPhaseDiffPerOcc(lowSpeedIdxes));
             
              meanPhaseOffsetHighSpeed = circMeanDeg(avgPhasePerOcc(highSpeedIdxes));
               meanPhaseOffsetLowSpeed = circMeanDeg(avgPhasePerOcc(lowSpeedIdxes));
            
               if(nonCircStats)
                   masterPhaseDiffHighMinusLowSpeed=[masterPhaseDiffHighMinusLowSpeed;(meanPhaseDiffHighSpeed- meanPhaseDiffLowSpeed) ];
                 masterMeanPhaseHighMinusLowSpeed=[masterMeanPhaseHighMinusLowSpeed; (meanPhaseOffsetHighSpeed-meanPhaseOffsetLowSpeed) ];

               else
                 masterPhaseDiffHighMinusLowSpeed=[masterPhaseDiffHighMinusLowSpeed; angdiffDeg([meanPhaseDiffLowSpeed meanPhaseDiffHighSpeed])];
                 masterMeanPhaseHighMinusLowSpeed=[masterMeanPhaseHighMinusLowSpeed; angdiffDeg([meanPhaseOffsetLowSpeed meanPhaseOffsetHighSpeed])];
               end
            
            
            masterMeanDistFracHighMinusLowSpeed=[masterMeanDistFracHighMinusLowSpeed; (meanDistFracHighSpeed-meanDistFracLowSpeed)];
            %{
            if(showPlots)
                 figure; 
                 numImages=length(currSeqData.currThetaSeqSummaryImgPaths);
                for ai=1:numImages
                    subplot(floor(sqrt(numImages)),floor(sqrt(numImages))+1,ai)
                   imshow(currSeqData.currThetaSeqSummaryImgPaths{ai})
                end
                maxFig
                close all
            end
            %}
        end
        
        
    end
    
end
%%
fH=figure
phaseDiffEdges=linspace(-180, 180, 71);
circMeasureSpreadEdges=linspace(0,0.2,71);

if(zScoreSpeed)
    %travSpeedEdges=linspace(-3,3,21);
    travSpeedEdges=linspace(minZspeed,maxZspeed,21);
 
else
    travSpeedEdges=linspace(minSpeed,maxSpeed,21);
end
phaseEdges=linspace(0,360,71);

withSmooth=1;

subplot(1,2,1)
%{
if(useCircMeasureSpread)
    masterPhaseDiff(masterPhaseDiff<0.0001)=NaN;
    getJointDistrGivenX(masterTravSpeed,masterPhaseDiff,travSpeedEdges,circMeasureSpreadEdges,fH,withSmooth)
else
    getJointDistrGivenX(masterTravSpeed,masterPhaseDiff,travSpeedEdges,phaseDiffEdges,fH,withSmooth)

end
%}
getJointDistrGivenX(masterTravSpeed,masterPhaseDiff,travSpeedEdges,phaseDiffEdges,fH,withSmooth)
%getJointDistrGivenXpeakNorm(masterTravSpeed,masterPhaseDiff,travSpeedEdges,phaseDiffEdges,fH,withSmooth)

xlabel('Traversal speed (m/s)')
 if(speedInFieldsPerSec) 
     xlabel('Traversal speed (fields/sec)')
 end
ylabel('Avg phase diff in theta sequence (deg)')
axis square
title(sprintf('Speed vs phase separation, n=%d theta sequences',length(masterMeanPhase)))


subplot(1,2,2)
getJointDistrGivenX(masterTravSpeed,masterMeanPhase,travSpeedEdges,phaseEdges,fH,withSmooth)
%getJointDistrGivenXpeakNorm(masterTravSpeed,masterMeanPhase,travSpeedEdges,phaseEdges,fH,withSmooth)

axis square

xlabel('Traversal speed (m/s)')
 if(speedInFieldsPerSec) 
     xlabel('Traversal speed (fields/sec)')
 end
ylabel('Avg phase in theta sequence (deg)')
title(sprintf('Speed vs phase offset, n=%d theta sequences',length(masterMeanPhase)))
setFigFontTo(16)
maxFig
saveas(gcf,'SpeedVsThetaSequencePhaseDiffAndOffset.png')
%%
%getBinnedCircAverages
if(zScoreSpeed)
  %travSpeedAvgEdges=linspace(-3,3,16);  
  travSpeedAvgEdges=linspace(minZspeed,maxZspeed,16);  
else
travSpeedAvgEdges=linspace(minSpeed,maxSpeed,16);
end
%travSpeedAvgEdges=linspace(minSpeed,maxSpeed,15);
%[speedBinCenters,phaseAvgValues] = getBinnedCircAverages(masterTravSpeed,masterMeanPhase,travSpeedAvgEdges);
%[speedBinCenters,phaseDiffValues] = getBinnedCircAverages(masterTravSpeed,masterPhaseDiff,travSpeedAvgEdges);
[speedBinCenters,phaseAvgValues] = getBinnedCircAverages(masterTravSpeed,masterMeanPhase,travSpeedAvgEdges);

travSpeedAvgDiffEdges=linspace(minSpeed,maxSpeed,16);

[speedBinDiffCenters,phaseDiffValues] = getBinnedCircAverages(masterTravSpeed,masterPhaseDiff,travSpeedAvgDiffEdges);



figure;
subplot(1,2,1)
plot(speedBinCenters,phaseAvgValues,'k.','MarkerSize',30)
xlabel('Traversal speed (m/s)')

if(speedInFieldsPerSec)
    xlabel('Traversal speed (fields/sec)')
end
ylabel('Avg phase in theta sequence (deg)')

ylim([140 240])
if(zScoreSpeed)
    xlim([minZspeed maxZspeed])
else
    xlim([minSpeed maxSpeed])
end
axis square
box off

subplot(1,2,2)
plot(speedBinDiffCenters,phaseDiffValues,'k.','MarkerSize',30)
xlabel('Traversal speed (m/s)')
if(speedInFieldsPerSec)
    xlabel('Traversal speed (fields/sec)')
end
ylabel('Avg phase diff in theta sequence (deg)')
ylim([0 100])
if(zScoreSpeed)
    xlim([minZspeed maxZspeed])
else
    xlim([minSpeed maxSpeed])
end
axis square
box off

saveas(gcf,fullfile(offsetAndDiffVsSpeedDir,'SpeedVsThetaSequencePhaseDiffAndOffset_BinnedAvgs.png'))
%figure; histogram(masterPhaseDiffSpeedRhos,100)
%figure; histogram(masterMeanPhaseSpeedRho,100)

%%

if(isempty(masterPhaseDiffHighMinusLowSpeed))
    
end 

[pMeanDiff,h,statsMeanDiff] = signrank(masterPhaseDiffHighMinusLowSpeed);

[pOffset,h,statsOffset] = signrank(masterMeanPhaseHighMinusLowSpeed);

gaussKernelWidthDeg=10;
gaussKernelWidthDeg=15;
%gaussKernelWidthDeg=5;
figure; 
phaseDiffEdges=linspace(-180,180,101);
phaseDiffBins=edgesToBins(phaseDiffEdges);
subplot(1,2,2);
%histogram(masterPhaseDiffHighMinusLowSpeed,phaseDiffEdges)

[pDistPhaseDiff,~]=getCircKernelDistr(masterPhaseDiffHighMinusLowSpeed,phaseDiffEdges,gaussKernelWidthDeg);
plot(phaseDiffBins,pDistPhaseDiff,'k-','LineWidth',4)
xlim([-180 180])
ylim([0 max(pDistPhaseDiff)])
%ylim([0 0.027])
hold on
xlabel('high speed circ-minus low speed phase difference (deg)')
ylabel('Probability')
vline(0,'k--',3)
title({'Speed vs phase separation of individual field sequence representations',sprintf('Wilcoxon signed rank test, z=%.2f, p=%.5f',statsMeanDiff.zval,pMeanDiff)})

box off
subplot(1,2,1);
%histogram(masterMeanPhaseHighMinusLowSpeed,phaseDiffEdges)
[pDistPhaseOffset,~]=getCircKernelDistr(masterMeanPhaseHighMinusLowSpeed,phaseDiffEdges,gaussKernelWidthDeg);
%plot(phaseDiffBins,pDistPhaseOffset,'k-','LineWidth',4)
plot(phaseDiffBins,pDistPhaseOffset,'r-','LineWidth',4)
xlim([-180 180])

ylim([0 max(pDistPhaseOffset)])


xlabel('high speed circ-minus low speed mean phase (deg)')
ylabel('Probability')
hold on
vline(0,'k--',3)
vline(circMeanDegNoMod360(masterMeanPhaseHighMinusLowSpeed(:)),'r--',3)
title({'Speed vs phase offset of individual field sequence representations',sprintf('Wilcoxon signed rank test, z=%.2f, p=%.5f',statsOffset.zval,pOffset)})
setFigFontTo(18)
maxFigHalfWidth
box off
maxFigMukkaalWidth
saveas(gcf,fullfile(offsetAndDiffVsSpeedDir,'indThetaSeqPhaseOffsetVsPhaseDiffsStats.png'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%no significant difference in field fraction between high and low speed gps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; histogram(masterMeanDistFracHighMinusLowSpeed,21)
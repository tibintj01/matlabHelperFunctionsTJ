close all;  clc
%dataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles';

clearvars -except spikeDataPerField

%offsetAndDiffVsSpeedDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3Images/offsetAndDiffVsSpeedIndThetaSeqMin2FieldsCntrlDist';
offsetAndDiffVsSpeedDir='individualThetaSequencesAnyTrajLengthMinNumFields_4_WithFiringRate';

%nonCircStats=1;
nonCircStats=0;

controlForDistanceAcrossLaps=1;

touchDir(offsetAndDiffVsSpeedDir)
tic
disp('loading spike data per field')
if(~exist('spikeDataPerField','var'))
    spikeDataPerField=load('allUnitInFieldGoodSpeedSpikePhasesAndTimes.mat');
end
toc

maxRunSpeed=1;

processedDataDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData';
unitDataDir=fullfile(processedDataDir,'unitSpikeInfoStructs');
dataDir=unitDataDir;

%saveThetaSeqImgDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3Images/individualThetaSequencesSimilarSizeMinTrajLengthHalfWidth';
saveThetaSeqImgDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3Images/individualThetaSequencesAnyTrajLengthMinNumFields_4_WithFiringRate';

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

%minSpeedStdDevAcrossLaps=0.1;
minSpeedStdDevAcrossLaps=0;
%maxSpeedStdDevWithinTrav=0.1;
maxSpeedStdDevWithinTrav=1;

minNumPtsInSpeedCategory=2;
minNumPtsInSpeedCategory=1;

%minSpeedStdDevAcrossLaps=0.15;
minNumOccurences=5;
minNumFieldsInSeq=2;
minNumFieldsInSeq=4;
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

indFieldHighSpeedPrctileCutoff=80;
indFieldLowSpeedPrctileCutoff=20;

maxAllowableSlopeMagnitude=720;
maxAllowableSlopeMagnitude=Inf;

minSpeed=0.1;
minSpeed=0.125;
%maxSpeed=0.7;
maxSpeed=0.7;

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
    
    if(~strcmp(currSessName,'ec013.626'))
        continue
    end
    
    %{
    if(isempty(currThetaSeqData.thetaSequenceInfoPerSequence))
        continue
    end
    %}
    
    if(~isfield(currThetaSeqData,'thetaSequenceInfoPerSequenceAtLeast4Fields'))
        continue
    end
    
    %thetaSequenceInfoPerSequence=currThetaSeqData.thetaSequenceInfoPerSequence;
        thetaSequenceInfoPerSequence=currThetaSeqData.thetaSequenceInfoPerSequenceAtLeast4Fields;

    
    seqStrs=fieldnames(thetaSequenceInfoPerSequence);
    
    for ti=1:length(seqStrs)
        currSeqData=thetaSequenceInfoPerSequence.(seqStrs{ti});
        numOccurences=currSeqData.numOccurences;
        numFieldsInSeq=length(currSeqData.fieldCenterSeq);
        
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
        
        
        
            peakPhaseSeqPerOcc=currSeqData.peakPhaseSeqPerOcc;
            travAvgSpeedPerOcc=currSeqData.travAvgSpeedPerOcc;
            travSpeedTracePerOcc=currSeqData.travSpeedTracePerOcc;
            numFieldsInSeq=length(currSeqData.fieldCenterSeq);
            
            avgPhaseDiffPerOcc=NaN(numOccurences,1);
            avgPhasePerOcc=NaN(numOccurences,1);
            speedStdDevWithinTrav=NaN(numOccurences,1);
            
            for oi=1:numOccurences
                
                if(~isControlledDistPerOcc(oi))
                    continue
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
            travSpeedPerOcc=cell2num(travAvgSpeedPerOcc);
            
            computeCorrIdxes=travSpeedPerOcc>=minSpeed & travSpeedPerOcc<=maxSpeed;
            
            constSpeedTravIdxes=speedStdDevWithinTrav <= maxSpeedStdDevWithinTrav;
            
            avgPhaseDiffPerOcc(~constSpeedTravIdxes)=NaN;
            avgPhasePerOcc(~constSpeedTravIdxes)=NaN;
            
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
            
            %highSpeedIdxes=computeCorrIdxes & travSpeedPerOcc>=highSpeedMinCutoff;
            %lowSpeedIdxes=computeCorrIdxes & travSpeedPerOcc<=lowSpeedMaxCutoff;
            
             %highSpeedIdxes=travSpeedPerOcc>=highSpeedMinCutoff;
            %lowSpeedIdxes=travSpeedPerOcc<=lowSpeedMaxCutoff;
            
             highSpeedIdxes=constSpeedTravIdxes(:) & travSpeedPerOcc(:)>=highSpeedMinCutoff;
            lowSpeedIdxes=constSpeedTravIdxes(:) & travSpeedPerOcc(:)<=lowSpeedMaxCutoff;
            
            %figure place holder
            circFigH=[];
            
            enoughPtsInEachCategory=sum(highSpeedIdxes)>=minNumPtsInSpeedCategory && sum(lowSpeedIdxes)>=minNumPtsInSpeedCategory;
            if(showPlots && enoughPtsInEachCategory)
                subplot(1,2,1)
                travSpeedPerOcc=cell2num(travAvgSpeedPerOcc);
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
            
            if(showPlots && enoughPtsInEachCategory)
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
           
            
            [currSeqPhaseDiffSpeedRho,currSeqPhaseDiffSpeedPval,currSeqPhaseDiffSpeedSlopeDegPerMsec, ~]...
                =getCircCorrCoeff(travSpeedPerOcc(computeCorrIdxes),avgPhaseDiffPerOcc(computeCorrIdxes),circFigH);
            
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
            [currSeqMeanPhaseSpeedRho,currSeqMeanPhaseSpeedPval,currSeqMeanPhaseSpeedSlopeDegPerMsec, ~]...
                =getCircCorrCoeff(travSpeedPerOcc(computeCorrIdxes),avgPhasePerOcc(computeCorrIdxes),circFigH);
            
            
               
            
            
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
            
           
                  
             if(abs(currSeqPhaseDiffSpeedSlopeDegPerMsec)<= maxAllowableSlopeMagnitude ...
                     && abs(currSeqMeanPhaseSpeedSlopeDegPerMsec)<=maxAllowableSlopeMagnitude)
                 
                 masterPhaseDiffSpeedRhos=[masterPhaseDiffSpeedRhos; currSeqPhaseDiffSpeedRho(:)];
                 masterMeanPhaseSpeedRho=[masterMeanPhaseSpeedRho; currSeqMeanPhaseSpeedRho(:)];
          
                 masterPhaseDiffSpeedSlope=[masterPhaseDiffSpeedSlope; currSeqPhaseDiffSpeedSlopeDegPerMsec(:)];
                 masterMeanPhaseSpeedSlope=[masterMeanPhaseSpeedSlope; currSeqMeanPhaseSpeedSlopeDegPerMsec(:)];
             

                masterTravSpeed=[masterTravSpeed; travSpeedPerOcc(:)];
                masterPhaseDiff=[masterPhaseDiff; avgPhaseDiffPerOcc(:)];
                masterMeanPhase=[masterMeanPhase; avgPhasePerOcc(:)];
            
             end
             
             if(sum(highSpeedIdxes)<minNumPtsInSpeedCategory || sum(lowSpeedIdxes)<minNumPtsInSpeedCategory)
                continue
             end
            
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %get high speed and low speed values
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
travSpeedEdges=linspace(minSpeed,maxSpeed,21);
phaseEdges=linspace(0,360,71);
withSmooth=1;

subplot(1,2,1)
if(useCircMeasureSpread)
    masterPhaseDiff(masterPhaseDiff<0.0001)=NaN;
    getJointDistrGivenX(masterTravSpeed,masterPhaseDiff,travSpeedEdges,circMeasureSpreadEdges,fH,withSmooth)
else
    getJointDistrGivenX(masterTravSpeed,masterPhaseDiff,travSpeedEdges,phaseDiffEdges,fH,withSmooth)

end
xlabel('Traversal speed (m/s)')
ylabel('Avg phase diff in theta sequence (deg)')
axis square
title(sprintf('Speed vs phase separation, n=%d theta sequences',length(masterMeanPhase)))


subplot(1,2,2)
getJointDistrGivenX(masterTravSpeed,masterMeanPhase,travSpeedEdges,phaseEdges,fH,withSmooth)
axis square

xlabel('Traversal speed (m/s)')
ylabel('Avg phase in theta sequence (deg)')
title(sprintf('Speed vs phase offset, n=%d theta sequences',length(masterMeanPhase)))
setFigFontTo(16)
maxFig
saveas(gcf,'SpeedVsThetaSequencePhaseDiffAndOffset.png')
%%
%getBinnedCircAverages
travSpeedAvgEdges=linspace(minSpeed,maxSpeed,16);
%travSpeedAvgEdges=linspace(minSpeed,maxSpeed,15);
[speedBinCenters,phaseAvgValues] = getBinnedCircAverages(masterTravSpeed,masterMeanPhase,travSpeedAvgEdges);
[speedBinCenters,phaseDiffValues] = getBinnedCircAverages(masterTravSpeed,masterPhaseDiff,travSpeedAvgEdges);



figure;
subplot(1,2,1)
plot(speedBinCenters,phaseAvgValues,'k.','MarkerSize',30)
xlabel('Traversal speed (m/s)')
ylabel('Avg phase in theta sequence (deg)')

ylim([140 240])
xlim([minSpeed maxSpeed])
axis square
box off

subplot(1,2,2)
plot(speedBinCenters,phaseDiffValues,'k.','MarkerSize',30)
xlabel('Traversal speed (m/s)')
ylabel('Avg phase diff in theta sequence (deg)')
ylim([0 100])
xlim([minSpeed maxSpeed])
axis square
box off

%figure; histogram(masterPhaseDiffSpeedRhos,100)
%figure; histogram(masterMeanPhaseSpeedRho,100)

%%

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
ylim([0 0.03])
hold on
xlabel('high speed circ-minus low speed phase difference (deg)')
ylabel('Probability')
vline(0,'k--',3)
title({'Speed vs phase separation of individual field sequence representations',sprintf('Wilcoxon signed rank test, z=%.2f, p=%.4f',statsMeanDiff.zval,pMeanDiff)})

box off
subplot(1,2,1);
%histogram(masterMeanPhaseHighMinusLowSpeed,phaseDiffEdges)
[pDistPhaseOffset,~]=getCircKernelDistr(masterMeanPhaseHighMinusLowSpeed,phaseDiffEdges,gaussKernelWidthDeg);
%plot(phaseDiffBins,pDistPhaseOffset,'k-','LineWidth',4)
plot(phaseDiffBins,pDistPhaseOffset,'r-','LineWidth',4)
xlim([-180 180])

ylim([0 0.03])
xlabel('high speed circ-minus low speed mean phase (deg)')
ylabel('Probability')
hold on
vline(0,'k--',3)
vline(circMeanDegNoMod360(masterMeanPhaseHighMinusLowSpeed(:)),'r--',3)
title({'Speed vs phase offset of individual field sequence representations',sprintf('Wilcoxon signed rank test, z=%.2f, p=%.4f',statsOffset.zval,pOffset)})
setFigFontTo(18)
maxFigHalfWidth
box off
saveas(gcf,fullfile(offsetAndDiffVsSpeedDir,'indThetaSeqPhaseOffsetVsPhaseDiffsStats.png'))
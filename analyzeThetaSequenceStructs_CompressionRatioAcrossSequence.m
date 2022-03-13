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
    
    
    
    if(isempty(currThetaSeqData.thetaSequenceInfoPerSequence))
        continue
    end
    
    
    thetaSequenceInfoPerSequence=currThetaSeqData.thetaSequenceInfoPerSequence;


    
    seqStrs=fieldnames(thetaSequenceInfoPerSequence);
    
    for ti=1:length(seqStrs)
        currSeqData=thetaSequenceInfoPerSequence.(seqStrs{ti});
        numOccurences=currSeqData.numOccurences;
       
        avgFieldWidth=currSeqData.avgFieldWidth;
        numFieldsInSeq=length(currSeqData.fieldCenterSeq);
               
        seqDirStr=currSeqData.seqDirStr;
            
        allPhaseSeqPerOriginalOcc=currSeqData.allPhaseSeqPerOcc;
        distAlongTrackPerOcc=currSeqData.positionPerOcc;
        fieldCenterDistAlongTrackPerField=currSeqData.fieldCenterSeq;
        
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
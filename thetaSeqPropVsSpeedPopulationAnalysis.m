close all;  clc
clearvars -except spikeDataPerField
offsetAndDiffVsSpeedDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3Images/offsetAndDiffVsSpeedIndThetaSeqMin2FieldsCntrlDist';

phaseSepColor=[255 208 0 ]/255;
phaseTimeCompColor=[109 197 165]/255;
meanPhaseColor=[0 77 138]/255;

plotPeakNorm=1;
plotPeakNorm=0;
zScoreSpeed=0;
thetaSeqPropVsSpeedInfo=load('thetaSequencePropVsSpeedProcessedData.mat');

masterTravSpeed=thetaSeqPropVsSpeedInfo.allSessionRunSpeed;
masterMeanDistFrac=thetaSeqPropVsSpeedInfo.allSessionMeanDistFracs;
masterDistSepFrac=thetaSeqPropVsSpeedInfo.allSessionAvgDistSepFracs;

masterPhaseDiff=thetaSeqPropVsSpeedInfo.allSessionPhaseSep;
masterMeanPhase=thetaSeqPropVsSpeedInfo.allSessionMeanPhase;
masterPhaseTimeComps=thetaSeqPropVsSpeedInfo.allSessionPhaseComp;

lowSpeedMax=0.55; %field/sec
highSpeedMin=0.95; %field/sec

lowSpeedMax=0.65; %field/sec
highSpeedMin=1; %field/sec

lowSpeedMax=0.6; %field/sec
highSpeedMin=1.1; %field/sec

lowSpeedMax=0.55; %field/sec
highSpeedMin=1.05; %field/sec


%lowSpeedMax=0.6; %field/sec
%highSpeedMin=1.2; %field/sec

maxMeanDistFrac=0.95;
minMeanDistFrac=0.05;


maxMeanDistFrac=0.7;
minMeanDistFrac=0.3;

maxMeanDistFrac=0.65;
minMeanDistFrac=0.35;


maxMeanDistFrac=0.75;
minMeanDistFrac=0.25;


maxMeanDistFrac=0.9;
minMeanDistFrac=0.1;

%{
maxMeanDistFrac=1;
minMeanDistFrac=0;
%}

minLowSpeedMeanPhase=0;

%{
maxMeanDistFrac=0.8;
minMeanDistFrac=0.2;

maxMeanDistFrac=0.7;
minMeanDistFrac=0.3;
%}



%maxMeanDistFrac=1;
%minMeanDistFrac=0;
%{
maxMeanDistFrac=0.85;
minMeanDistFrac=0.15;
%}


%{
maxDistRangeFrac=0.25;
maxDistRangeFrac=0.3;
maxDistRangeFrac=0.2;
%maxDistRangeFrac=0.25;
maxDistRangeFrac=0.3;
%}

maxDistRangeFrac=0.4;
maxDistRangeFrac=0.45;
%maxDistRangeFrac=0.6;
%maxDistRangeFrac=0.35;
minDistRangeFrac=0.1;
%minDistRangeFrac=0.05;


maxDistRangeFrac=0.4;
minDistRangeFrac=0.1;

maxDistRangeFrac=0.7;
minDistRangeFrac=0.1;

maxDistRangeFrac=0.85;
minDistRangeFrac=0.15;
%minDistRangeFrac=0.2;



%{
maxDistRangeFrac=1;
minDistRangeFrac=0;
%}


%maxDistRangeFrac=0.5;
%minDistRangeFrac=0.;


validIdxes=masterMeanDistFrac<=maxMeanDistFrac & masterMeanDistFrac>=minMeanDistFrac & masterDistSepFrac<=maxDistRangeFrac & masterDistSepFrac>=minDistRangeFrac;

masterTravSpeed(~validIdxes)=NaN;
masterMeanDistFrac(~validIdxes)=NaN;
masterDistSepFrac(~validIdxes)=NaN;

masterPhaseDiff(~validIdxes)=NaN;
masterMeanPhase(~validIdxes)=NaN;
masterPhaseTimeComps(~validIdxes)=NaN;

maxCompression=1080;
maxCompression=720; %deg per second
%maxCompression=480; %deg per second

%maxCompression=480;

masterPhaseTimeComps(masterPhaseTimeComps>maxCompression)=NaN;


   minSpeed=0.2;
   maxSpeed=1.5;
   
   
   minSpeed=0.25;
      maxSpeed=1.4;
   
        minSpeed=0.4;
      maxSpeed=1.2;
      
       minSpeed=0.4;
      maxSpeed=1.1;
      
         minSpeed=0.25;
      maxSpeed=1.4;
      
       minSpeed=0.3;
      maxSpeed=1.2;
      
       minSpeed=0.2;
   maxSpeed=1.2;
   
    minSpeed=0.5;
   maxSpeed=1.25;
   
    minSpeed=0.4;
   maxSpeed=1.2;
   
     minSpeed=0.25;
   maxSpeed=1.25;
   
       minSpeed=0.4;
   maxSpeed=1.25;
   
     minSpeed=0.2;
   maxSpeed=1.4;
   
   minSpeed=0.25;
   minSpeed=0.3;
   maxSpeed=1.4;
   maxSpeed=1.35;
   
   
   
    minSpeed=0.15;
    %minSpeed=0.1;
     %minSpeed=0.2;
   %maxSpeed=1.2;
   maxSpeed=1.45;
   
    minSpeed=0.2;
   maxSpeed=1.2;
   
   minSpeed=0.25;
   maxSpeed=1.25;
   
    minSpeed=0.25;
   maxSpeed=1.25;
   
   minSpeed=0.25;
   maxSpeed=1.2;
   
      minSpeed=0.25;
   maxSpeed=1.25;
   
    minSpeed=0.3;
   maxSpeed=1.2;
   
    minSpeed=0.3;
   maxSpeed=1.25;
 

      %minSpeed=0.2;
      

   %maxSpeed=1.2;
   
   %minSpeed=0.15;
   %maxSpeed=1.5;
   %maxSpeed=1.2;
   
   %{
   minSpeed=0.4;
   maxSpeed=1.3;
   %}
   
fH=figure
phaseDiffEdges=linspace(-180, 180, 31);

travSpeedEdges=linspace(minSpeed,maxSpeed,21);

travSpeedEdges=linspace(minSpeed,maxSpeed,41);

phaseEdges=linspace(0,360,51);
phaseTimeCompEdges=linspace(0,540,31);

phaseTimeCompEdges=linspace(0,720,101);
distFracEdges=linspace(0,1,51);

numRunningAvgBins=10;
numRunningAvgBins=8;
numRunningAvgBins=15;
numRunningAvgBins=12;
%numRunningAvgBins=20;
%numRunningAvgBins=15;
%numRunningAvgBins=12;

%numRunningAvgBins=15;
%numRunningAvgBins=12;
numRunningAvgBins=16;
%numRunningAvgBins=20;
numRunningAvgBins=10;
numRunningAvgBins=20;
numRunningAvgBins=16;
numRunningAvgBins=12;
numRunningAvgBins=20;
numRunningAvgBins=16;
numRunningAvgBins=14;
numRunningAvgBins=10;
numRunningAvgBins=14;
numRunningAvgBins=12;
numRunningAvgBins=20;
numRunningAvgBins=10;
numRunningAvgBins=12;

travSpeedAvgEdges=linspace(minSpeed,maxSpeed,numRunningAvgBins+1);

%travSpeedAvgEdges=linspace(minSpeed,maxSpeed,15);
%[speedBinCenters,phaseAvgValues] = getBinnedCircAverages(masterTravSpeed,masterMeanPhase,travSpeedAvgEdges);
%[speedBinCenters,phaseDiffValues] = getBinnedCircAverages(masterTravSpeed,masterPhaseDiff,travSpeedAvgEdges);

masterLogTravSpeed=getLogTime(masterTravSpeed);

masterTravSpeedOriginal=masterTravSpeed;
%masterTravSpeed=masterLogTravSpeed;
[speedBinCenters,phaseAvgValues,lowerLimMeanPhase,upperLimMeanPhase] = getBinnedCircAverages(masterTravSpeed,masterMeanPhase,travSpeedAvgEdges);

travSpeedAvgDiffEdges=linspace(minSpeed,maxSpeed,numRunningAvgBins+1);

[speedBinDiffCenters,phaseDiffValues,lowerLimPhaseSep,upperLimPhaseSep] = getBinnedCircAverages(masterTravSpeed,masterPhaseDiff,travSpeedAvgDiffEdges);

[runningAvgSpeedBinCentersComp,distFracRunningAvgValues,distSEMperSpeedBin] = getBinnedAverages(masterTravSpeed,masterMeanDistFrac,travSpeedAvgEdges);

[runningAvgSpeedBinCentersComp,distFracRangeRunningAvgValues,distRangeSEMperSpeedBin] = getBinnedAverages(masterTravSpeed,masterDistSepFrac,travSpeedAvgEdges);


phaseCompValidIdxes=masterPhaseTimeComps<=maxCompression & masterPhaseTimeComps>=0 & ~isnan(masterPhaseTimeComps);
[runningAvgSpeedBinCentersComp,phaseTimeCompRunningAvgValues,phaseTimeCompSEMperSpeedBin,phaseTimeCompStdDevPerSpeedBin] = getBinnedAverages(masterTravSpeed(phaseCompValidIdxes),masterPhaseTimeComps(phaseCompValidIdxes),travSpeedAvgEdges);

withSmooth=1;

fD=figure;
subplot(2,2,1)
getJointDistrGivenX(masterTravSpeed,masterMeanDistFrac,travSpeedEdges,distFracEdges,fD,withSmooth)

subplot(2,2,2)
plot(runningAvgSpeedBinCentersComp,distFracRunningAvgValues,'k.','MarkerSize',30)
hold on
shadedErrorBar(runningAvgSpeedBinCentersComp,distFracRunningAvgValues,distSEMperSpeedBin)
xlabel('Traversal speed (fields/s)')
ylabel('Avg distance in field (frac)')

subplot(2,2,3)
getJointDistrGivenX(masterTravSpeed,masterDistSepFrac,travSpeedEdges,distFracEdges,fD,withSmooth)

subplot(2,2,4)
plot(runningAvgSpeedBinCentersComp,distFracRangeRunningAvgValues,'k.','MarkerSize',30)
hold on
shadedErrorBar(runningAvgSpeedBinCentersComp,distFracRangeRunningAvgValues,distRangeSEMperSpeedBin)

%%
lowSpeedMasterIdxes=masterTravSpeed<=lowSpeedMax;
highSpeedMasterIdxes=masterTravSpeed>=highSpeedMin;

meanDistLowSpeedMean=nanmean(masterMeanDistFrac(lowSpeedMasterIdxes));
meanDistHighSpeedMean=nanmean(masterMeanDistFrac(highSpeedMasterIdxes));

meanDistLowSpeedSEM=nanstd(masterMeanDistFrac(lowSpeedMasterIdxes));
meanDistHighSpeedSEM=nanstd(masterMeanDistFrac(highSpeedMasterIdxes));

rangeDistLowSpeedMean=nanmean(masterDistSepFrac(lowSpeedMasterIdxes));
rangeDistHighSpeedMean=nanmean(masterDistSepFrac(highSpeedMasterIdxes));

rangeDistLowSpeedSEM=nanstd(masterDistSepFrac(lowSpeedMasterIdxes));
rangeDistHighSpeedSEM=nanstd(masterDistSepFrac(highSpeedMasterIdxes));

figure;

subplot(1,2,1)
bar([meanDistLowSpeedMean meanDistHighSpeedMean],'w')
hold on
errorbar(1,meanDistLowSpeedMean,meanDistLowSpeedSEM,'k')
errorbar(2,meanDistHighSpeedMean,meanDistHighSpeedSEM,'k')

title('Distance in field')
ylabel('Mean distance in field (frac)')
%xlim([0.5 1.5])
box off

subplot(1,2,2)
bar([rangeDistLowSpeedMean rangeDistHighSpeedMean],'w')
hold on
errorbar(1,rangeDistLowSpeedMean,rangeDistLowSpeedSEM,'k')
errorbar(2,rangeDistHighSpeedMean,rangeDistHighSpeedSEM,'k')

title('Distance separation between fields')
ylabel('Distance separation (frac)')

box off
%xlim([0.5 1.5])

setFigFontTo(18)
maxFigHalfHalfWidth

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%distance control stats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%figure; histogram(masterMeanDistFrac(lowSpeedMasterIdxes)); yyaxis right; histogram(masterMeanDistFrac(highSpeedMasterIdxes))

fH=figure

subplot(1,3,1)

if(plotPeakNorm)
    getJointDistrGivenXpeakNorm(masterTravSpeed,masterPhaseDiff,travSpeedEdges,phaseDiffEdges,fH,withSmooth)

else
    getJointDistrGivenX(masterTravSpeed,masterPhaseDiff,travSpeedEdges,phaseDiffEdges,fH,withSmooth)
end

     xlabel('Traversal speed (fields/sec)')
 
ylabel('Avg phase diff in theta sequence (deg)')
axis square
title(sprintf('Speed vs phase separation, n=%d theta sequences',sum(validIdxes)))



subplot(1,3,2)
if(plotPeakNorm)
    getJointDistrGivenXpeakNorm(masterTravSpeed,masterMeanPhase,travSpeedEdges,phaseEdges,fH,withSmooth)
else
    getJointDistrGivenX(masterTravSpeed,masterMeanPhase,travSpeedEdges,phaseEdges,fH,withSmooth)
end

axis square

     xlabel('Traversal speed (fields/sec)')

ylabel('Avg phase in theta sequence (deg)')
title(sprintf('Speed vs phase offset, n=%d theta sequences',sum(validIdxes)))

%withSmooth=1;

subplot(1,3,3)
if(plotPeakNorm)
    getJointDistrGivenXpeakNorm(masterTravSpeed,masterPhaseTimeComps,travSpeedEdges,phaseTimeCompEdges,fH,withSmooth)
else
    getJointDistrGivenX(masterTravSpeed,masterPhaseTimeComps,travSpeedEdges,phaseTimeCompEdges,fH,withSmooth)
end

axis square

xlabel('Traversal speed (fields/sec)')

ylabel('Phase-time compression ratio (deg/field duration)')
title(sprintf('Speed vs phase-time compression, n=%d theta sequences',sum(validIdxes)))



setFigFontTo(16)
maxFig
%saveas(gcf,'SpeedVsThetaSequencePhaseDiffAndOffset.png')

if(plotPeakNorm)
    saveas(gcf,fullfile(offsetAndDiffVsSpeedDir,'SpeedVsThetaSequencePhaseDiffAndOffset_PeakNormDistr.png'))
else
    saveas(gcf,fullfile(offsetAndDiffVsSpeedDir,'SpeedVsThetaSequencePhaseDiffAndOffset.png'))
end



%%
figure;
subplot(1,3,3)
plot(speedBinCenters,phaseAvgValues,'.','MarkerSize',30,'Color',meanPhaseColor)
hold on
%lowerLimMeanPhase
negErr=(phaseAvgValues-lowerLimMeanPhase);
posErr=(upperLimMeanPhase-phaseAvgValues);

errorbar(speedBinCenters,phaseAvgValues, negErr, posErr,'.', 'Color',meanPhaseColor)


[m,b,R,p]=getLinearFit(runningAvgSpeedBinCentersComp,phaseAvgValues,minSpeed,maxSpeed,meanPhaseColor);
title(sprintf('Linear regression R=%.2f, p=%.2d',R,p))
    xlabel('Traversal speed (fields/sec)')

ylabel('Avg phase in theta sequence (deg)')

ylim([75 275])
ylim([120 270])
ylim([140 290])

ylim([120 300])
ylim([130 290])

%ylim([150 270])
if(zScoreSpeed)
    xlim([minZspeed maxZspeed])
else
    xlim([minSpeed maxSpeed])
end
axis square
box off

subplot(1,3,1)
plot(speedBinDiffCenters,phaseDiffValues,'k.','MarkerSize',30,'Color',phaseSepColor)
hold on
negErr=(phaseDiffValues-lowerLimPhaseSep);
posErr=(upperLimPhaseSep-phaseDiffValues);

errorbar(speedBinDiffCenters,phaseDiffValues, negErr, posErr,'.','Color',phaseSepColor)
[m,b,R,p]=getLinearFit(runningAvgSpeedBinCentersComp,phaseDiffValues,minSpeed,maxSpeed,phaseSepColor);
title(sprintf('Linear regression R=%.2f, p=%.2f',R,p))

    xlabel('Traversal speed (fields/sec)')

ylabel('Avg phase diff in theta sequence (deg)')
%ylim([-100 100])
ylim([-30 120]+20)

ylim([0 180]-20)
%ylim([150 270]-140)
if(zScoreSpeed)
    xlim([minZspeed maxZspeed])
else
    xlim([minSpeed maxSpeed])
end
axis square
box off

subplot(1,3,2)
plot(runningAvgSpeedBinCentersComp,phaseTimeCompRunningAvgValues,'.','MarkerSize',30,'Color',phaseTimeCompColor)
hold on
errorbar(runningAvgSpeedBinCentersComp,phaseTimeCompRunningAvgValues, phaseTimeCompSEMperSpeedBin,'.','Color',phaseTimeCompColor)
hold on
[m,b,R,p]=getLinearFit(runningAvgSpeedBinCentersComp,phaseTimeCompRunningAvgValues,minSpeed,maxSpeed,phaseTimeCompColor);

title(sprintf('Linear regression R=%.2f, p=%.2d',R,p))

    xlabel('Traversal speed (fields/sec)')

ylabel('Phase-time compression ratio (deg/sec)')
%ylim([-100 100])
ylim([180 480])
ylim([180 400])
ylim([200 380])
ylim([180 330])
ylim([150 340])

    xlim([minSpeed maxSpeed])

axis square
box off

%phaseTimeCompSEMperSpeedBin
maxFig
setFigFontTo(16)

masterTravSpeed=masterTravSpeedOriginal;

saveas(gcf,fullfile(offsetAndDiffVsSpeedDir,'SpeedVsThetaSequencePhaseDiffAndOffset_BinnedAvgs.png'))
%figure; histogram(masterPhaseDiffSpeedRhos,100)
%figure; histogram(masterMeanPhaseSpeedRho,100)

%%

distSepPerFieldSequence=thetaSeqPropVsSpeedInfo.distSepPerFieldSequence;
timeSepPerFieldSequence=thetaSeqPropVsSpeedInfo.timeSepPerFieldSequence;
meanDistPerFieldSequence=thetaSeqPropVsSpeedInfo.meanDistPerFieldSequence;
runSpeedsPerFieldSequence=thetaSeqPropVsSpeedInfo.runSpeedsPerFieldSequence;
meanPhasePerFieldSequence=thetaSeqPropVsSpeedInfo.meanPhasePerFieldSequence;
phaseSepPerFieldSequence=thetaSeqPropVsSpeedInfo.phaseSepPerFieldSequence;
phaseCompPerFieldSequence=thetaSeqPropVsSpeedInfo.phaseCompPerFieldSequence;
fieldSeqCount =thetaSeqPropVsSpeedInfo.fieldSeqCount;

masterPhaseDiffHighMinusLowSpeed=[];
masterMeanPhaseHighMinusLowSpeed=[];
masterPhaseDiffHighMinusLowSpeed=[];
masterPhaseCompHighMinusLowSpeed=[];

masterPhaseDiffHighSpeed=[];
         masterPhaseDiffLowSpeed=[];
         
         masterMeanPhaseHighSpeed=[];
         masterMeanPhaseLowSpeed=[];
         masterMeanSpeedDifferential=[];
         
         masterPhaseCompHighSpeed=[];
         masterPhaseCompLowSpeed=[];
         
         masterMeanDistHighSpeed=[];
         masterMeanDistLowSpeed=[];
         
         masterDistSepHighSpeed=[];
         masterDistSepLowSpeed=[];
         
           masterTimeSepHighSpeed=[];
         masterTimeSepLowSpeed=[];
         

         figure
for fsi=1:fieldSeqCount
    currSeqRunSpeeds=runSpeedsPerFieldSequence{fsi};
    currSeqMeanPhases=meanPhasePerFieldSequence{fsi};
    currSeqPhaseSeps=phaseSepPerFieldSequence{fsi};
    currSeqPhaseComps=phaseCompPerFieldSequence{fsi};
    
    currSeqMeanDists=meanDistPerFieldSequence{fsi};
    currSeqDistSep=distSepPerFieldSequence{fsi};
    currSeqTimeSep=timeSepPerFieldSequence{fsi};
    
    %{
    validIdxes=currSeqMeanDists<=maxMeanDistFrac & currSeqMeanDists>=minMeanDistFrac & currSeqDistSep<=maxDistRangeFrac & currSeqDistSep>=minDistRangeFrac;

    
    currSeqRunSpeeds(~validIdxes)=NaN;
    currSeqMeanDists(~validIdxes)=NaN;
    currSeqDistSep(~validIdxes)=NaN;

    currSeqPhaseSeps(~validIdxes)=NaN;
    currSeqMeanPhases(~validIdxes)=NaN;
    currSeqPhaseComps(~validIdxes)=NaN;
    %}
    
    
    if(nanmean(currSeqDistSep)>maxDistRangeFrac || nanmean(currSeqDistSep) <minDistRangeFrac ...
            ||  nanmean(currSeqMeanDists)>maxMeanDistFrac  || nanmean(currSeqMeanDists)<minMeanDistFrac)
            continue
            
    end
    
    
    
    %{
    plot([prctile(currSeqRunSpeeds,25) prctile(currSeqRunSpeeds,75)],[fsi fsi],'k-')
    hold on
    continue
    %}
    
    currSeqLowSpeedIdxes=currSeqRunSpeeds<=lowSpeedMax & currSeqRunSpeeds>=0;
    currSeqHighSpeedIdxes=currSeqRunSpeeds>=highSpeedMin & currSeqRunSpeeds<=Inf;
    
     %currSeqLowSpeedIdxes=currSeqRunSpeeds<=lowSpeedMax & currSeqRunSpeeds>=minSpeed;
    %currSeqHighSpeedIdxes=currSeqRunSpeeds>=highSpeedMin & currSeqRunSpeeds<=maxSpeed;
    

    %{
    currSeqLowSpeedMax=prctile(currSeqRunSpeeds,25);
    currSeqHighSpeedMin=prctile(currSeqRunSpeeds,75);
    
    currSeqLowSpeedIdxes=currSeqRunSpeeds<=currSeqLowSpeedMax & currSeqRunSpeeds>=minSpeed;
    currSeqHighSpeedIdxes=currSeqRunSpeeds>=currSeqHighSpeedMin & currSeqRunSpeeds<=maxSpeed;
    %}
    
    minNumPtsInCategory=1;
     %minNumPtsInCategory=3;
    if(sum(currSeqLowSpeedIdxes)>=minNumPtsInCategory && sum(currSeqHighSpeedIdxes)>=minNumPtsInCategory)
        %highSpeedPhaseSep=circMeanDegNoMod360(currSeqPhaseSeps(currSeqHighSpeedIdxes));
        %lowSpeedPhaseSep=circMeanDegNoMod360(currSeqPhaseSeps(currSeqLowSpeedIdxes));
        
        highSpeedPhaseSep=nanmean(currSeqPhaseSeps(currSeqHighSpeedIdxes));
        lowSpeedPhaseSep=nanmean(currSeqPhaseSeps(currSeqLowSpeedIdxes));

        %highSpeedMeanPhase=circMeanDeg(currSeqMeanPhases(currSeqHighSpeedIdxes));
        %lowSpeedMeanPhase=circMeanDeg(currSeqMeanPhases(currSeqLowSpeedIdxes));
        
           highSpeedMeanPhase=nanmean(currSeqMeanPhases(currSeqHighSpeedIdxes));
        lowSpeedMeanPhase=nanmean(currSeqMeanPhases(currSeqLowSpeedIdxes));
        
        highSpeedMeanDist=nanmean(currSeqMeanDists(currSeqHighSpeedIdxes));
        lowSpeedMeanDist=nanmean(currSeqMeanDists(currSeqLowSpeedIdxes));
        
        currSeqMeanSpeedDiff=nanmean(currSeqRunSpeeds(currSeqHighSpeedIdxes))-nanmean(currSeqRunSpeeds(currSeqLowSpeedIdxes));

        highSpeedDistSep=nanmean(currSeqDistSep(currSeqHighSpeedIdxes));
        lowSpeedDistSep=nanmean(currSeqDistSep(currSeqLowSpeedIdxes));
        
        highSpeedTimeSep=nanmean(currSeqTimeSep(currSeqHighSpeedIdxes));
        lowSpeedTimeSep=nanmean(currSeqTimeSep(currSeqLowSpeedIdxes));
        
        
        
        
        highSpeedPhaseComp=abs(nanmean(currSeqPhaseComps(currSeqHighSpeedIdxes)));
        lowSpeedPhaseComp=abs(nanmean(currSeqPhaseComps(currSeqLowSpeedIdxes)));
        
        if(lowSpeedMeanPhase>minLowSpeedMeanPhase)

            %masterIndSeqSpeeds=[masterIndSeqSpeeds; 
             masterPhaseDiffHighSpeed=[masterPhaseDiffHighSpeed; highSpeedPhaseSep];
             masterPhaseDiffLowSpeed=[masterPhaseDiffLowSpeed; lowSpeedPhaseSep];

             masterMeanPhaseHighSpeed=[masterMeanPhaseHighSpeed; highSpeedMeanPhase];
             masterMeanPhaseLowSpeed=[masterMeanPhaseLowSpeed; lowSpeedMeanPhase];
             masterMeanSpeedDifferential=[masterMeanSpeedDifferential; currSeqMeanSpeedDiff];

             masterPhaseCompHighSpeed=[masterPhaseCompHighSpeed; highSpeedPhaseComp];
             masterPhaseCompLowSpeed=[masterPhaseCompLowSpeed; lowSpeedPhaseComp];

              masterMeanDistHighSpeed=[masterMeanDistHighSpeed; highSpeedMeanDist];
             masterMeanDistLowSpeed=[masterMeanDistLowSpeed; lowSpeedMeanDist];

             masterDistSepHighSpeed=[masterDistSepHighSpeed; highSpeedDistSep];
             masterDistSepLowSpeed=[masterDistSepLowSpeed; lowSpeedDistSep];
             
              masterTimeSepHighSpeed=[masterTimeSepHighSpeed; highSpeedTimeSep];
             masterTimeSepLowSpeed=[masterTimeSepLowSpeed; lowSpeedTimeSep];
             
             

            masterPhaseDiffHighMinusLowSpeed=[masterPhaseDiffHighMinusLowSpeed; angdiffDeg( [lowSpeedPhaseSep highSpeedPhaseSep ])];
            masterMeanPhaseHighMinusLowSpeed=[masterMeanPhaseHighMinusLowSpeed; angdiffDeg([lowSpeedMeanPhase highSpeedMeanPhase ])];
            masterPhaseCompHighMinusLowSpeed=[masterPhaseCompHighMinusLowSpeed; (highSpeedPhaseComp- lowSpeedPhaseComp)];
        
        end
    end
    
end
%%
%distance mean and separation stats
Nbootstrap=100000;
%Nbootstrap=10000;
%[surrogateMeanDiffs,originalMeanDiff]=getBootstrapTwoGpCircMeanDiff(masterMeanPhaseLowSpeed,masterMeanPhaseHighSpeed,Nbootstrap);
[surrogateMeanDistDiffs,originalMeanDistDiff]=getBootstrapTwoGpMeanDiff(masterMeanDistLowSpeed,masterMeanDistHighSpeed,Nbootstrap);
%[surrogateMeanDiffs,originalMeanDiff]=getBootstrapTwoGpMeanCircDiff(masterMeanPhaseLowSpeed,masterMeanPhaseHighSpeed,Nbootstrap);

figure; %phaseDiffEdges=linspace(-180,180,100);
subplot(2,3,1)
meanDistDiffEdges=linspace(-0.1,0.1,200);
%[currPdist,~]=getCircKernelDistr(surrogateMeanDiffs,phaseDiffEdges,1);
%[currPdist,~]=getKernelDistrAndPeakVal(surrogateMeanDiffs,0.5,phaseDiffEdges);
%getProbDist(vals,edges,showPlots,isCirc)
[currPdist]=getProbDist(surrogateMeanDistDiffs,meanDistDiffEdges,0,0);

meanDistDiffBins=edgesToBins(meanDistDiffEdges);
plot(meanDistDiffBins,currPdist,'k-')
hold on
axis tight
vline(originalMeanDistDiff)
box off

[pVal]=getBootstrapPVal(currPdist,meanDistDiffBins,originalMeanDistDiff);

 xlabel('Mean distance (high speed - low speed) (frac)'); ylabel('Probability')
 title(sprintf('Bootstrap shuffled p value: %.4f',pVal))
 
 
Nbootstrap=100000;
%[surrogateDiffDiffs,originalDiffDiff]=getBootstrapTwoGpCircMeanDiff(masterPhaseDiffLowSpeed,masterPhaseDiffHighSpeed,Nbootstrap);
[surrogateDistDiffDiffs,originalDistDiffDiff]=getBootstrapTwoGpMeanDiff(masterDistSepLowSpeed,masterDistSepHighSpeed,Nbootstrap);
%[surrogateDiffDiffs,originalDiffDiff]=getBootstrapTwoGpMeanCircDiff(masterPhaseDiffLowSpeed,masterPhaseDiffHighSpeed,Nbootstrap);


[surrogateTimeDiffDiffs,originalTimeDiffDiff]=getBootstrapTwoGpMeanDiff(masterTimeSepLowSpeed,masterTimeSepHighSpeed,Nbootstrap);

subplot(2,3,2)
%phaseDiffEdges=linspace(-180,180,100);
distSepDiffEdges=linspace(-0.1,0.1,200);
%[currPdist,~]=getCircKernelDistr(surrogateMeanDiffs,phaseDiffEdges,1);
%[currPdist,~]=getKernelDistrAndPeakVal(surrogateDiffDiffs,1,phaseDiffEdges);
%getProbDist(vals,edges,showPlots,isCirc)
[currPdist]=getProbDist(surrogateDistDiffDiffs,distSepDiffEdges,0,0);

distSepDiffBins=edgesToBins(distSepDiffEdges);
plot(distSepDiffBins,currPdist,'k-')
hold on
axis tight
vline(originalDistDiffDiff)
box off

[pVal]=getBootstrapPVal(currPdist,distSepDiffBins,originalDistDiffDiff);

 xlabel('Avg distance separation (high speed - low speed) (frac)'); ylabel('Probability')
 title(sprintf('Bootstrap shuffled p value: %.4f',pVal))
 
 
  subplot(2,3,3)
 
 
 %phaseDiffEdges=linspace(-180,180,100);
timeSepDiffEdges=linspace(-0.5,0.5,200);
%[currPdist,~]=getCircKernelDistr(surrogateMeanDiffs,phaseDiffEdges,1);
%[currPdist,~]=getKernelDistrAndPeakVal(surrogateDiffDiffs,1,phaseDiffEdges);
%getProbDist(vals,edges,showPlots,isCirc)
[currPdist]=getProbDist(surrogateTimeDiffDiffs,timeSepDiffEdges,0,0);

timeSepDiffBins=edgesToBins(timeSepDiffEdges);
plot(timeSepDiffBins,currPdist,'k-')
hold on
axis tight
box off
vline(originalTimeDiffDiff)

[pVal]=getBootstrapPVal(currPdist,timeSepDiffBins,originalTimeDiffDiff);

 xlabel('Avg time separation (high speed - low speed) (frac)'); ylabel('Probability')
 title(sprintf('Bootstrap shuffled p value: %d',pVal))
 
 %bar plots
 subplot(2,3,4)
 
 meanDistLowSpeedMean=nanmean(masterMeanDistLowSpeed);
 meanDistLowSpeedSEM=getSEMacrossRows(masterMeanDistLowSpeed(:));
 meanDistHighSpeedMean=nanmean(masterMeanDistHighSpeed);
 meanDistHighSpeedSEM=getSEMacrossRows(masterMeanDistHighSpeed(:));
 
bar([meanDistLowSpeedMean meanDistHighSpeedMean],'w')
hold on
errorbar(1,meanDistLowSpeedMean,meanDistLowSpeedSEM,'k')
errorbar(2,meanDistHighSpeedMean,meanDistHighSpeedSEM,'k')

title('Distance in field')
ylabel('Mean distance in field (frac)')
%xlim([0.5 1.5])
box off

 rangeDistLowSpeedMean=nanmean(masterDistSepLowSpeed);
 rangeDistLowSpeedSEM=getSEMacrossRows(masterDistSepLowSpeed(:));
 rangeDistHighSpeedMean=nanmean(masterDistSepHighSpeed);
 rangeDistHighSpeedSEM=getSEMacrossRows(masterDistSepHighSpeed(:));

 subplot(2,3,5)
bar([rangeDistLowSpeedMean rangeDistHighSpeedMean],'w')
hold on
errorbar(1,rangeDistLowSpeedMean,rangeDistLowSpeedSEM,'k')
errorbar(2,rangeDistHighSpeedMean,rangeDistHighSpeedSEM,'k')

title('Distance separation between fields')
ylabel('Distance separation (frac)')

box off
 

 rangeTimeLowSpeedMean=nanmean(masterTimeSepLowSpeed);
 rangeTimeLowSpeedSEM=getSEMacrossRows(masterTimeSepLowSpeed(:));
 rangeTimeHighSpeedMean=nanmean(masterTimeSepHighSpeed);
 rangeTimeHighSpeedSEM=getSEMacrossRows(masterTimeSepHighSpeed(:));

 subplot(2,3,6)
bar([rangeTimeLowSpeedMean rangeTimeHighSpeedMean],'w')
hold on
errorbar(1,rangeTimeLowSpeedMean,rangeTimeLowSpeedSEM,'k')
errorbar(2,rangeTimeHighSpeedMean,rangeTimeHighSpeedSEM,'k')

title('Time separation between fields')
ylabel('Time separation (frac)')

box off
 
  

 
  setFigFontTo(18)
  maxFig
 saveas(gcf,'distanceAndTimeBootstrapStats.png')

%%
figure;
%colormap(jet)

subplot(1,3,1)

ptSize=200;
scatter(masterPhaseDiffLowSpeed,masterPhaseDiffHighSpeed,ptSize,phaseSepColor,'filled'); axis square
%colorbar
hold on; plot([-180 180], [-180 180],'k--','LineWidth',3)
xlim([-180 180])
ylim([-180 180])
xlabel('Low speed phase separation (deg)')
ylabel('High speed phase separation (deg)')

maxCompDisplay=500;
%maxCompDisplay=720;
subplot(1,3,2)
scatter(masterPhaseCompLowSpeed,masterPhaseCompHighSpeed,ptSize,phaseTimeCompColor,'filled'); axis square
%colorbar
hold on; plot([0 maxCompDisplay], [0 maxCompDisplay],'k--','LineWidth',3)
xlim([0 maxCompDisplay])
ylim([0 maxCompDisplay])
xlabel('Low speed phase-time ratio (deg/sec)')
ylabel('High speed phase-time ratio (deg/sec)')

subplot(1,3,3)
scatter(masterMeanPhaseLowSpeed,masterMeanPhaseHighSpeed,ptSize,meanPhaseColor,'filled'); axis square
%colorbar
hold on; plot([0 360], [0 360],'k--','LineWidth',3)
xlim([0 360])
ylim([0 360])
xlabel('Low speed mean phase (deg)')
ylabel('High speed mean phase (deg)')

setFigFontTo(18)
maxFig
%caxis([0.6 0.9])
saveas(gcf,'lowSpeedToHighSpeedRelativeToXEqualsY.png')

%{
highOriginalPhaseIdxes=masterMeanPhaseLowSpeed>180;
masterMeanPhaseLowSpeed(~highOriginalPhaseIdxes)=NaN;
masterMeanPhaseHighSpeed(~highOriginalPhaseIdxes)=NaN;
%}
%%
Nbootstrap=500000;
%[surrogateMeanDiffs,originalMeanDiff]=getBootstrapTwoGpCircMeanDiff(masterMeanPhaseLowSpeed,masterMeanPhaseHighSpeed,Nbootstrap);
[surrogateMeanDiffs,originalMeanDiff]=getBootstrapTwoGpMeanDiff(masterMeanPhaseLowSpeed,masterMeanPhaseHighSpeed,Nbootstrap);
%[surrogateMeanDiffs,originalMeanDiff]=getBootstrapTwoGpMeanCircDiff(masterMeanPhaseLowSpeed,masterMeanPhaseHighSpeed,Nbootstrap);

figure; %phaseDiffEdges=linspace(-180,180,100);
subplot(1,3,3)
phaseDiffEdges=linspace(-50,50,200);
%[currPdist,~]=getCircKernelDistr(surrogateMeanDiffs,phaseDiffEdges,1);
%[currPdist,~]=getKernelDistrAndPeakVal(surrogateMeanDiffs,0.5,phaseDiffEdges);
%getProbDist(vals,edges,showPlots,isCirc)
[currPdist]=getProbDist(surrogateMeanDiffs,phaseDiffEdges,0,0);


lineWidthBootstrap=2;
originalMeanLineWidth=5;
phaseDiffBins=edgesToBins(phaseDiffEdges);
plot(phaseDiffBins,currPdist,'-','Color','k','LineWidth',lineWidthBootstrap)
hold on
axis tight
vline(originalMeanDiff,meanPhaseColor,originalMeanLineWidth)


[pVal]=getBootstrapPVal(currPdist,phaseDiffBins,originalMeanDiff);

 xlabel('Mean phase (high speed - low speed) (deg)'); ylabel('Probability')
 title(sprintf('Bootstrap shuffled p value: %.4f',pVal))
 
 %setFigFontTo(18)
 %saveas(gcf,'meanPhaseBootstrapStats.png')
 

%[surrogateDiffDiffs,originalDiffDiff]=getBootstrapTwoGpCircMeanDiff(masterPhaseDiffLowSpeed,masterPhaseDiffHighSpeed,Nbootstrap);
[surrogateDiffDiffs,originalDiffDiff]=getBootstrapTwoGpMeanDiff(masterPhaseDiffLowSpeed,masterPhaseDiffHighSpeed,Nbootstrap);
%[surrogateDiffDiffs,originalDiffDiff]=getBootstrapTwoGpMeanCircDiff(masterPhaseDiffLowSpeed,masterPhaseDiffHighSpeed,Nbootstrap);

 box off

subplot(1,3,1)
%phaseDiffEdges=linspace(-180,180,100);
phaseDiffEdges=linspace(-50,50,200);
%[currPdist,~]=getCircKernelDistr(surrogateMeanDiffs,phaseDiffEdges,1);
%[currPdist,~]=getKernelDistrAndPeakVal(surrogateDiffDiffs,1,phaseDiffEdges);
%getProbDist(vals,edges,showPlots,isCirc)
[currPdist]=getProbDist(surrogateDiffDiffs,phaseDiffEdges,0,0);


phaseDiffBins=edgesToBins(phaseDiffEdges);
plot(phaseDiffBins,currPdist,'-','Color','k','LineWidth',lineWidthBootstrap)
hold on
axis tight
vline(originalDiffDiff,phaseSepColor,originalMeanLineWidth)


[pVal]=getBootstrapPVal(currPdist,phaseDiffBins,originalDiffDiff);

 xlabel('Avg phase separation (high speed - low speed) (deg)'); ylabel('Probability')
 title(sprintf('Bootstrap shuffled p value: %.4f',pVal))
 

%[surrogateDiffDiffs,originalDiffDiff]=getBootstrapTwoGpCircMeanDiff(masterPhaseDiffLowSpeed,masterPhaseDiffHighSpeed,Nbootstrap);
%validCompIdxes=masterPhaseCompLowSpeed<=maxCompression & masterPhaseCompHighSpeed<=maxCompression ;
validCompIdxes=masterPhaseCompLowSpeed<=Inf & masterPhaseCompHighSpeed<=Inf ;


[surrogateCompDiffs,originalCompDiff]=getBootstrapTwoGpMeanDiff(masterPhaseCompLowSpeed(validCompIdxes),masterPhaseCompHighSpeed(validCompIdxes),Nbootstrap);
 box off 

 subplot(1,3,2)
%phaseDiffEdges=linspace(-180,180,100);
%phaseCompEdges=linspace(-500,500,200);
phaseCompEdges=linspace(-400,400,200);
%[currPdist,~]=getCircKernelDistr(surrogateMeanDiffs,phaseDiffEdges,1);
%[currPdist,~]=getKernelDistrAndPeakVal(surrogateCompDiffs,10,phaseCompEdges);
%getProbDist(vals,edges,showPlots,isCirc)
[currPdist]=getProbDist(surrogateCompDiffs,phaseCompEdges,0,0);


phaseCompBins=edgesToBins(phaseCompEdges);
plot(phaseCompBins,currPdist,'-','Color','k','LineWidth',lineWidthBootstrap)
hold on
axis tight
vline(originalCompDiff,phaseTimeCompColor,originalMeanLineWidth)


[pVal]=getBootstrapPVal(currPdist,phaseCompBins,originalCompDiff);

 xlabel('Avg phase-time comp (high speed - low speed) (deg)'); ylabel('Probability')
 title(sprintf('Bootstrap shuffled p value: %.4f',pVal))
 
 box off
 
  setFigFontTo(18)
  maxFig
 saveas(gcf,'speedPhaseSigBootstrapStats.png')

%%
%{
[surrogateDiffDiffs,originalDiffDiff]=getBootstrapTwoGpCircMeanDiff(masterPhaseDiffLowSpeed,masterPhaseDiffHighSpeed,Nbootstrap);

figure; phaseDiffEdges=linspace(-180,180,100);[currPdist,~]=getCircKernelDistr(surrogateDiffDiffs,phaseDiffEdges,1); plot(edgesToBins(phaseDiffEdges),currPdist,'k-')

[pMeanDistDiff,h,statsMeanDistDiff] = signrank(masterDistSepHighSpeed,masterDistSepLowSpeed)

[pMeanDist,h,statsMeanDist] = signrank(masterMeanDistHighSpeed,masterMeanDistLowSpeed)

%}


if(isempty(masterPhaseDiffHighMinusLowSpeed))
    
end 

[pMeanDiff,h,statsMeanDiff] = signrank(masterPhaseDiffHighMinusLowSpeed);

[pOffset,h,statsOffset] = signrank(masterMeanPhaseHighMinusLowSpeed);

[pComp,h,statsComp] = signrank(masterPhaseCompHighMinusLowSpeed);

gaussKernelWidthDeg=10;
gaussKernelWidthDeg=15;
gaussKernelWidthDeg=20;
%gaussKernelWidthDeg=5;
figure; 
phaseDiffEdges=linspace(-180,180,101);
phaseDiffBins=edgesToBins(phaseDiffEdges);
subplot(1,3,1);
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
title({'Speed vs phase separation per field sequence',sprintf('Wilcoxon signed rank test, z=%.2f, p=%.5f',statsMeanDiff.zval,pMeanDiff)})

box off
subplot(1,3,2);
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
title({'Speed vs phase offset per field sequence',sprintf('Wilcoxon signed rank test, z=%.2f, p=%.5f',statsOffset.zval,pOffset)})
setFigFontTo(18)
maxFigHalfWidth
box off

subplot(1,3,3);
maxCompDiff=1000;
phaseTimeCompDiffEdges=linspace(-maxCompDiff, maxCompDiff, 31);
phaseTimeCompDiffBins=edgesToBins(phaseTimeCompDiffEdges);
%histogram(masterMeanPhaseHighMinusLowSpeed,phaseDiffEdges)

[pDistPhaseComp,~]=getKernelDistrAndPeakVal(masterPhaseCompHighMinusLowSpeed,100,phaseTimeCompDiffEdges);
%plot(phaseDiffBins,pDistPhaseOffset,'k-','LineWidth',4)
plot(phaseTimeCompDiffBins,pDistPhaseComp,'r-','LineWidth',4)
xlim([-maxCompDiff maxCompDiff])

ylim([0 max(pDistPhaseComp)])


xlabel('high speed minus low speed phase-time ratio (deg/sec)')
ylabel('Probability')
hold on
vline(0,'k--',3)
validIdxes=abs(masterPhaseCompHighMinusLowSpeed)<=maxCompDiff;

vline(nanmean(masterPhaseCompHighMinusLowSpeed(validIdxes)),'r--',3)
title({'Speed vs phase comp per field sequence',sprintf('Wilcoxon signed rank test, z=%.2f, p=%.5f',statsComp.zval,pComp)})
setFigFontTo(16)
maxFig
box off

%maxFigMukkaalWidth
saveas(gcf,fullfile(offsetAndDiffVsSpeedDir,'indThetaSeqPhaseOffsetVsPhaseDiffsStats.png'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%no significant difference in field fraction between high and low speed gps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure; histogram(masterMeanDistFracHighMinusLowSpeed,21)
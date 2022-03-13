%another property of log compression is linear relationship between time and width 
close all; clear all; clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load time and phase data within good field traversals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
setTightSubplots_Spacious
data=load('allUnitInFieldGoodSpeedSpikePhasesAndTimes.mat')

saveFigureDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/figurePanels';
touchDir(saveFigureDir)

%totalUsedFieldCount=data.totalUsedFieldCount;
totalFieldCount=data.totalFieldCount;

masterSpikePhases=data.allUnitSpikePhases;
%masterFieldTimes=data.allUnitSpikeTimesInField/0.8;
masterFieldTimes=data.allUnitSpikeTimesInField/0.9;

timeOffset=0;
%masterFieldTimes=getLogTime((masterFieldTimes-timeOffset),1.333);

numTimeBins=10;
numTimeBins=8;
numTimeBins=7;
numTimeBins=12;
numTimeBins=50;
numTimeBins=7;
numTimeBins=10;
numTimeBins=8;
numTimeBins=12;
timeEdges=linspace(0,1,numTimeBins+1);
timeBinCenters=edgesToBins(timeEdges);
timeBinWidth=median(diff(timeBinCenters));


%numPhaseBins=24;
numPhaseBins=30;
numPhaseBins=60;
%numPhaseBins=90;
phaseEdges=linspace(0,360,numPhaseBins+1);
phaseBinCenters=edgesToBins(phaseEdges);
base=1.333;
%masterFieldTimes=getLogTime(1-masterFieldTimes,base);

jfH=figure;
[pxySmooth,pxy]=getJointDistr(masterFieldTimes,masterSpikePhases,timeEdges,phaseEdges,jfH);



%{
[~,phaseAvgValues] = getBinnedCircAverages(masterFieldTimes,masterSpikePhases,timeEdges);
figure; plot(timeBinCenters,phaseAvgValues,'ko-','LineWidth',4)
ylim([0 360])
%}
title(sprintf('n=%d fields',totalFieldCount))
xlabel('Time in field (frac)')
ylabel('Theta phase (deg)')


%daspect([1 360 1])
%caxis([0 prctile(pxySmooth(:),97.5)])
%maxFig
setFigFontTo(18)
saveas(gcf,'timeVsThetaPhaseJointDistr.tif')




timeBinColors=jet(numTimeBins)
maxProbPhasePerTimeBin=NaN(numTimeBins,1);
distrWidth=NaN(numTimeBins,1);
for ti=1:numTimeBins
    %currPdist=circSmooth(pxySmooth(ti,:),5);
    currPdist=circSmooth(pxy(ti,:),5);
    
    currPdist=currPdist/sum(currPdist);
    

      [maxProb,maxID]=max(currPdist);
%    currPdist=circshift(currPdist,-maxID);

subplot(1,2,1)
    
[phaseBinCentersPadded,currPdistPadded,padLength]=getCircPaddedDistr(phaseBinCenters,currPdist,10);

normFact=maxProb;
normFact=1;

    plot(phaseBinCentersPadded,currPdistPadded/normFact,'-','Color',timeBinColors(ti,:),'LineWidth',3)
    
   [currDistrWidth,risingIntersectPhase,risingIntersectProb,descendingIntersectPhase,descendingIntersectProb]...
       =getFullWidthHalfMaxCircDistr(phaseBinCenters,currPdist);
   
    distrWidth(ti)=currDistrWidth;
    hold on
    
    plot(risingIntersectPhase,risingIntersectProb/normFact,'.','Color',timeBinColors(ti,:),'MarkerSize',30)
    plot(descendingIntersectPhase,descendingIntersectProb/normFact,'.','Color',timeBinColors(ti,:),'MarkerSize',30)
    
    plot([risingIntersectPhase descendingIntersectPhase], [risingIntersectProb descendingIntersectProb]/normFact,'Color',timeBinColors(ti,:),'LineWidth',3)
  
    maxProbPhasePerTimeBin(ti)=phaseBinCenters(maxID);
    

    %circMeanPhasePerSpeedBin(ti)=circMeanDeg(
    
    
    hold on
    
    subplot(1,2,2)
    plot(timeBinCenters(ti),distrWidth(ti),'.','Color',timeBinColors(ti,:),'MarkerSize',80)
    hold on
end
subplot(1,2,1)
colormap(gca,jet)
cb=colorbar
ylabel(cb,'Time in field (frac)')
xlabel('Theta phase (deg)')

ylabel('Probability')
%xlim([0 360])
xlim([-20 380])
caxis([0 1])

title(sprintf('Theta phase distribution vs Time in traversal, n=%d fields',totalFieldCount))


    subplot(1,2,2)
    hold on
    %plot(timeBinCenters,distrWidth,'k-','LineWidth',3)
    
    ylim([100 340])
    
    xlabel('Time in field (frac)')
    ylabel('Theta phase distribution width')
    
    %rVal=getCorrCoeff(timeBinCenters,distrWidth)
    [m,b,R]=getLinearFit(timeBinCenters(1:(end)),distrWidth(1:(end)));
    
    title({sprintf('Distribution width vs time, R=%.2f',R),sprintf('slope=%.2f deg/field',m)})
    
    uberTitle({'Prediction of logarithmic time compression in theta cycles:', 'Linear increase in distribution width with time'})
    
    maxFig
setFigFontTo(18)

%saveas(gcf,'timeInFieldVsThetaPhaseDistr.tif')

saveToFigureDir('timeInFieldVsThetaPhaseDistrWidth',saveFigureDir)


close all; clear all;clc
tic
%load('speedVsOffsetData_Oct8_2021.mat')
load('speedVsOffsetData_HighResSpeed')
toc
close all


setTightSubplots
timeEdges=linspace(0,1,31);
phaseEdges=linspace(0,360,31);
allPhasesHiVsLowSpeedFig=figure;

offsetPerSpeed=NaN(numSpeedGroups,1);
slopePerSpeed=NaN(numSpeedGroups,1);
speedBinCentersDisp=NaN(numSpeedGroups,1);
base=1.3333;
base=1.4
%base=1.2
fH=figure;
timeOffset=0.2;


trueTimeBoundFact=0.8;
 

for si=1:numSpeedGroups
  currSpeedBinMin=speedCategoryBounds(si);
  currSpeedBinMax=speedCategoryBounds(si+1);
  speedBinCentersDisp(si)=(currSpeedBinMin+currSpeedBinMax)/2;
                      

  currSpeedGroupPhaseVsTime=allPhaseVsTimePerSpeedGroupAllFields{si};
  currGpTimeVals=currSpeedGroupPhaseVsTime(:,1)/trueTimeBoundFact;
  
  %currGpTimeVals=1-getLogTime(1-currSpeedGroupPhaseVsTime(:,1)-timeOffset,base);
  %currGpTimeVals=currGpTimeVals+normrnd(0.03,0.03,size(currGpTimeVals));
    

    
    if(true || si==1 || si==3 || si==5 || si==7)
     %   subplot(1,4,(si-1)/2+1)
        subplot(1,numSpeedGroups,si)
        getJointDistr(currGpTimeVals,currSpeedGroupPhaseVsTime(:,2),timeEdges,phaseEdges,allPhasesHiVsLowSpeedFig);
        %title(sprintf('speed < %.2f m/s',singleCellLowSpeedMax))
        title(sprintf('%.2f m/s < speed < %.2f m/s',currSpeedBinMin,currSpeedBinMax))
        [rho,p,slopeDegPerXunit, offsetDeg] = getCircCorrCoeff(currGpTimeVals,currSpeedGroupPhaseVsTime(:,2),fH)

        xlabel('Time in field (frac)')
        
        ylabel('Theta phase (deg)')
        caxis([0 2.5e-3])
    
    else
    
        [rho,p,slopeDegPerXunit, offsetDeg] = getCircCorrCoeff(currGpTimeVals,currSpeedGroupPhaseVsTime(:,2))
    end

    offsetPerSpeed(si)=offsetDeg;
    slopePerSpeed(si)=slopeDegPerXunit;
end
%%
figure; plot(speedBinCentersDisp, offsetPerSpeed,'k-o')

yyaxis right;
plot(speedBinCentersDisp, slopePerSpeed,'k-o')


%%

%close all; clear all;clc
tic
%load('speedVsOffsetData_CloseTimeBound.mat')
load('speedVsOffsetData_HighResSpeed.mat')
close all
toc

setTightSubplots_Medium

figure

subplot(1,3,1)
timeEdgesAvg=linspace(0,1,12);
%timeEdgesAvg=linspace(0,0.8,12);
%timeEdgesAvg=linspace(0.05,1,15);
%timeEdgesAvg=linspace(0,0.6,6);
%timeEdgesAvg=linspace(0,0.75,7);


speedGroupColors=jet(numSpeedGroups);

%speedGroupColors=redblue(numSpeedGroups);

compareGps=[3 7 ]
compareGps=[3 9 ]
compareGps=1:numSpeedGroups
%compareGps=1:9

localAvgOffsetPerSpeed=NaN(size(compareGps));
localAvgSlopePerSpeed=NaN(size(compareGps));
localAvgSpeedBinCenters=NaN(size(compareGps));

minPhaseDisp=70;


for si=compareGps
  currSpeedBinMin=speedCategoryBounds(si);
  currSpeedBinMax=speedCategoryBounds(si+1);
                      

  currSpeedGroupPhaseVsTime=allPhaseVsTimePerSpeedGroupAllFields{si};
  
  [timeBinCenters,currSpeedGpPhaseAvgValues] = getBinnedCircAverages(currSpeedGroupPhaseVsTime(:,1)/trueTimeBoundFact,currSpeedGroupPhaseVsTime(:,2),timeEdgesAvg);
  
  circShiftIdxesDisp=currSpeedGpPhaseAvgValues<minPhaseDisp;
  currSpeedGpPhaseAvgValues(circShiftIdxesDisp)=currSpeedGpPhaseAvgValues(circShiftIdxesDisp)+360;
 
  
  currSpeedGpPhaseAvgValues=smooth(currSpeedGpPhaseAvgValues,3);
    plot(timeBinCenters,currSpeedGpPhaseAvgValues,'-o','Color', speedGroupColors(si,:),'LineWidth',3)
    hold on
    
   [rho,p,slopeDegPerXunit, offsetDeg] = getCircCorrCoeff(timeBinCenters(:),currSpeedGpPhaseAvgValues(:));

            
 
    localAvgOffsetPerSpeed(si)=offsetDeg;
    localAvgSlopePerSpeed(si)=slopeDegPerXunit;
    localAvgSpeedBinCenters(si)=(currSpeedBinMax+currSpeedBinMin)/2;
end

colormap(jet)

cb=colorbar('north')

ylabel(cb,'Running speed (m/s)')
caxis([min(speedCategoryBounds) max(speedCategoryBounds)])
 xlabel('Time in field (frac)')
 ylabel('Circ mean theta phase (degrees)')
 box off
 setFigFontTo(16)
 
  title('Phase vs time vs speed')
subplot(1,3,2)
plot(localAvgSpeedBinCenters, localAvgOffsetPerSpeed,'k-','LineWidth',3)

%yyaxis right;
hold on
plot(localAvgSpeedBinCenters, localAvgOffsetPerSpeed,'k.','MarkerSize',50)


%plot(localAvgSpeedBinCenters, localAvgSlopePerSpeed,'k-o')
t=scaledata(linspace(0.35,0.7,100),0,0.85);

minOffset=-70;
maxOffset=-12;
minOffset=-75;
maxOffset=-22;

minOffset=-60;
maxOffset=-13;

minOffset=-75;
maxOffset=-28;

minOffset=-70;
maxOffset=-20;


minOffset=-80;
maxOffset=-5;

minOffset=-80;
maxOffset=-23;

minOffset=-84;
maxOffset=-28;

timeOffset=0.15;
logSpeed=scaledata(getLogTime(t+timeOffset,base),minOffset,maxOffset);


minSpeed=min(speedCategoryBounds);
maxSpeed=max(speedCategoryBounds);
speedRange=maxSpeed-minSpeed;
speedValues=scaledata(t,minSpeed,maxSpeed);

%speed in fields per second...

noRectify=1;
noLimit=1;
%%
figure;
testBases=1.1:0.1:3;

%testBases=exp(1);

for bi=1:length(testBases)
    base=testBases(bi);
    logSpeedBins=getLogTime((localAvgSpeedBinCenters-minSpeed)/speedRange+timeOffset,base,noRectify,noLimit);
    speedBinCenters
    subplot(4,5,bi)
     %plot(logSpeedBins,localAvgOffsetPerSpeed,'k.','MarkerSize',50)
     plot(localAvgSpeedBinCenters, localAvgOffsetPerSpeed,'k.','MarkerSize',50)
          %plot(t,logSpeed,'k.','MarkerSize',50)
          hold on
          logSpeed=scaledata(getLogTime(t+timeOffset,base),minOffset,maxOffset);
          plot(speedValues,logSpeed,'r--','LineWidth',4)


     
    [m,b,R,p]=getLinearFit(logSpeedBins,localAvgOffsetPerSpeed,0.2,1.1);
       % [m,b,R,p]=getLinearFit(logSpeedBins,localAvgOffsetPerSpeed);

    %xlim([0.2 1.1])
     xlim([0 1.1])
     %xlim([-Inf Inf])
    ylim([-90 -20])
    axis square

    xlabel(sprintf('log_{%.2f}(normalized speed)',base))
    ylabel('Phase precession offset (deg)')
    title(sprintf('Linear regression, R=%.2f, m=%.2f deg/(m/s), b=%.2f deg',R,m,b))
    box off

end
setFigFontTo(10)
%%

pL=plot(linspace(min(speedCategoryBounds),max(speedCategoryBounds),100),logSpeed,'r--','LineWidth',5)

xlabel('Running speed (m/s)')
ylabel('Phase precession offset (deg)')

 legend(pL,sprintf('log_{%.1f}(speed)',base),'Location','northwest')
 legend boxoff
 ylim([-100 -20])
 box off
 
 title('Phase-time offset vs speed')
 xlim([min(speedCategoryBounds),max(speedCategoryBounds)])
 subplot(1,3,3)
plot(localAvgSpeedBinCenters,localAvgSlopePerSpeed,'k.','MarkerSize',50)
hold on
plot(localAvgSpeedBinCenters,localAvgSlopePerSpeed,'k-','LineWidth',3)

hold on
plot(xlim,[nanmean(localAvgSlopePerSpeed) nanmean(localAvgSlopePerSpeed)],'r--','LineWidth',5)
title('Phase-time slope vs speed')
box off
xlim([min(speedCategoryBounds),max(speedCategoryBounds)])
 ylim([-270 0])
 xlabel('Running speed (m/s)')
ylabel('Phase precession slope (deg/field)')
  setFigFontTo(18)
  maxFig
  saveas(gcf,'runningSpeedVsOffsetLogBaseFit.png')

 %ylim([0 1])
 %%
autoArrangeFigures
%%
%collect 30cm/s vs 60cm/s slope and offset shuffled-resampling statistics
runShuffles=0;
if(runShuffles)
numShuffles=100000;

%lowSpeedCategoryIdx=2;
%highSpeedCategoryIdx=8;

lowSpeedCategoryIdx=2;
highSpeedCategoryIdx=9;

speedCategoryBounds(lowSpeedCategoryIdx)
speedCategoryBounds(highSpeedCategoryIdx)

%lowSpeedCategoryIdx=4;
%highSpeedCategoryIdx=6;

lowSpeedPhaseVsTime=allPhaseVsTimePerSpeedGroupAllFields{lowSpeedCategoryIdx};
highSpeedPhaseVsTime=allPhaseVsTimePerSpeedGroupAllFields{highSpeedCategoryIdx};


%{
lowSpeedOffset=offsetPerSpeed(lowSpeedCategoryIdx);
highSpeedOffset=offsetPerSpeed(highSpeedCategoryIdx);
lowSpeedSlope=slopePerSpeed(lowSpeedCategoryIdx);
highSpeedSlope=slopePerSpeed(highSpeedCategoryIdx);
%}

lowSpeedOffset=localAvgOffsetPerSpeed(lowSpeedCategoryIdx);
highSpeedOffset=localAvgOffsetPerSpeed(highSpeedCategoryIdx);
lowSpeedSlope=localAvgSlopePerSpeed(lowSpeedCategoryIdx);
highSpeedSlope=localAvgSlopePerSpeed(highSpeedCategoryIdx);

originalOffsetHighMinusLow=angdiffDeg([lowSpeedOffset highSpeedOffset])
%originalSlopeHighMinusLow=diff([lowSpeedSlope highSpeedSlope])
originalSlopeHighMinusLow=highSpeedSlope/lowSpeedSlope



combinedPhaseVsTime=[lowSpeedPhaseVsTime;highSpeedPhaseVsTime];

numLowSpeedPts=size(lowSpeedPhaseVsTime,1);
numHighSpeedPts=size(highSpeedPhaseVsTime,1);

shuffledOffsetDiffs=NaN(numShuffles,1);
shuffledSlopesRatios=NaN(numShuffles,1);

for shi=1:numShuffles
    if(mod(shi,100)==0)
        disp(shi)
    end
    
    lowSpeedSurrogatePtsIdxes=randsample(length(combinedPhaseVsTime),numLowSpeedPts);
    highSpeedSurrogatePtsIdxes=randsample(length(combinedPhaseVsTime),numHighSpeedPts);
    
    lowSpeedSurrogatePts=combinedPhaseVsTime(lowSpeedSurrogatePtsIdxes,:);
     highSpeedSurrogatePts=combinedPhaseVsTime(highSpeedSurrogatePtsIdxes,:);
     
    
     [timeBinCenters,lowSpeedGpPhaseAvgValues] = getBinnedCircAverages(lowSpeedSurrogatePts(:,1),lowSpeedSurrogatePts(:,2),timeEdgesAvg);
     [timeBinCenters,highSpeedGpPhaseAvgValues] = getBinnedCircAverages(highSpeedSurrogatePts(:,1),highSpeedSurrogatePts(:,2),timeEdgesAvg);

     
      [rho,p,lowSpeedSlopeDegPerSpeedChange, lowSpeedOffsetDeg] = getCircCorrCoeff(timeBinCenters(:),lowSpeedGpPhaseAvgValues(:));
     [rho,p,highSpeedSlopeDegPerSpeedChange, highSpeedOffsetDeg] = getCircCorrCoeff(timeBinCenters(:),highSpeedGpPhaseAvgValues(:));

     
    %{
    [rho,p,lowSpeedSlopeDegPerSpeedChange, lowSpeedOffsetDeg] = getCircCorrCoeff(lowSpeedSurrogatePts(:,1),lowSpeedSurrogatePts(:,2));
     [rho,p,highSpeedSlopeDegPerSpeedChange, highSpeedOffsetDeg] = getCircCorrCoeff(highSpeedSurrogatePts(:,1),highSpeedSurrogatePts(:,2));

     %}
     
     currShuffleOffsetHighMinusLow=angdiffDeg([lowSpeedOffsetDeg highSpeedOffsetDeg]);
	%currShuffleSlopeHighMinusLow=diff([lowSpeedSlopeDegPerSpeedChange highSpeedSlopeDegPerSpeedChange]);
    	currShuffleSlopeHighMinusLow=highSpeedSlopeDegPerSpeedChange/lowSpeedSlopeDegPerSpeedChange;

    
     
    shuffledOffsetDiffs(shi)=currShuffleOffsetHighMinusLow;
    shuffledSlopesRatios(shi)=currShuffleSlopeHighMinusLow;
     
end

%%
setTightSubplots_Spacious
close all
figure; subplot(2,1,1); 
%offsetDiffsEdges=linspace(-15,15,50);
offsetDiffsEdges=linspace(-15,35,100);
slopeRatiosEdges=linspace(0.8,1.25,50);

offsetDiffBins=edgesToBins(offsetDiffsEdges);
slopeRatioBins=edgesToBins(slopeRatiosEdges);

No=histcounts(shuffledOffsetDiffs,offsetDiffsEdges); 

plot(offsetDiffBins,No,'Color',[0.2 0 0],'LineWidth',4)
hold on;
ylim([0 max(No)])
plot([originalOffsetHighMinusLow originalOffsetHighMinusLow],ylim,'r--','LineWidth',5)
box off
xlabel('High speed to low speed offset difference (deg)')
ylabel('Shuffle count')
subplot(2,1,2); 
Ns=histcounts(shuffledSlopesDiffs,slopeRatiosEdges); 

plot(slopeRatioBins,Ns,'Color',[0.2 0 0],'LineWidth',4)
ylim([0 max(Ns)])
hold on;
plot([originalSlopeHighMinusLow originalSlopeHighMinusLow],ylim,'r--','LineWidth',5)
box off
xlabel('High speed to low speed slope ratio')
ylabel('Shuffle count')

setFigFontTo(18)
%save('speedVsOffsetStatsWithBootstrap.mat')
end





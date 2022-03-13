close all; clear all; clc
processedDataDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData';

saveDataPath=fullfile(processedDataDir,'phaseDistrWidthVsTimeInFieldData.mat');

data=load(saveDataPath);

allFieldTimeBinCenters=data.allFieldTimeBinCenters; %wrong name in previous code!

minR=-Inf;
%timeBinCenters=unique(allFieldTimeBinCenters);

%{
numPhaseBins=18;
phaseEdges=linspace(0,360,numPhaseBins+1);
phaseBinCenters=edgesToBins(phaseEdges);
%}

timeBinCenters=data.timeBinCenters;
allFieldDistrWidths=data.allFieldDistrWidths;

allFieldDistrWidths(allFieldDistrWidths>360)=NaN;

allFieldTimeWidthRs=data.allFieldTimeWidthRs;

numFieldsUsed=length(allFieldTimeWidthRs);

allFieldTimeWidthAvgs=data.allFieldTimeWidthAvgs;


allFieldRs=NaN(size(allFieldDistrWidths));
allFieldAvgWidths=NaN(size(allFieldDistrWidths));
for fi=1:numFieldsUsed
    currFieldValIdxes=(data.allFieldStatIDs==fi);
        
    allFieldRs(currFieldValIdxes)=allFieldTimeWidthRs(fi);
    allFieldAvgWidths(currFieldValIdxes)=allFieldTimeWidthAvgs(fi);
end

allValsHasGoodR=allFieldRs>minR;
allFieldDistrWidths(~allValsHasGoodR)=NaN;

allFieldLinearModelErrs=data.allFieldLinearModelErrs;
allFieldLinearModelErrs(~allValsHasGoodR)=NaN;


tolLev=0.01;

meanPhaseWidthsPerTimeBin=NaN(length(timeBinCenters),1);
meanPhaseWidthsErrPerTimeBin=NaN(length(timeBinCenters),1);

phaseErrSEM=NaN(length(timeBinCenters),1);

phaseSEM=NaN(length(timeBinCenters),1);


figure
for i=1:length(timeBinCenters)
    

    currTimeBinIdxes=(abs(allFieldTimeBinCenters-timeBinCenters(i))<tolLev);
    
    if(sum(currTimeBinIdxes)<=1)
        continue
    end
    
    currTimeBinLinModelErrs=allFieldLinearModelErrs(currTimeBinIdxes);
    %currTimeBinQuadModelErrs=allFieldQuadModelErrs(currTimeBinIdxes);
    
    currTimeBinWidths=allFieldDistrWidths(currTimeBinIdxes);
    
    currTimeBinWidthAvgs=allFieldAvgWidths(currTimeBinIdxes);
    %currTimeBinWidths=sqrt(currTimeBinWidths-120);
    %meanPhaseWidthsPerTimeBin(i)=circMeanDeg(currTimeBinWidths);
    %phaseSEM(i)=1/circ_kappa(currTimeBinWidths);
    
    %meanPhaseWidthsPerTimeBin(i)=nanmean(currTimeBinWidths-currTimeBinWidthAvgs);
        meanPhaseWidthsPerTimeBin(i)=nanmean(currTimeBinWidths);

    
    meanPhaseWidthsErrPerTimeBin(i)=nanmean(currTimeBinLinModelErrs);

    
    phaseErrSEM(i)=getSEMacrossRows(currTimeBinLinModelErrs(:));

    
    phaseSEM(i)=getSEMacrossRows(currTimeBinWidths(:));
    
        subplot(5,6,i)
        numWidthBins=50;
        widthEdges=linspace(0,360,numWidthBins+1);
        widthBins=edgesToBins(widthEdges);
        Nw=histcounts(currTimeBinWidths,widthEdges)
        %plot(widthBins,smooth(Nw,ceil(numWidthBins/5)),'k-','LineWidth',2)
        plot(widthBins,Nw,'k-','LineWidth',2)

        %xlim([100 200])
    
end

%meanPhaseWidthsPerTimeBin=sqrt(meanPhaseWidthsPerTimeBin);
setTightSubplots_Spacious

figure;
subplot(2,1,1)
plot(timeBinCenters,meanPhaseWidthsPerTimeBin,'ko')
hold on
shadedErrorBar(timeBinCenters,meanPhaseWidthsPerTimeBin,phaseSEM)
title('Theta phase distribution width vs time in field')

xlabel('Time in field (frac)')
%ylabel('Change in phase distr. width (deg)')
ylabel('Theta phase distribution width (deg)')
xlim([0 1])
ylim([65 165])

subplot(2,1,2)
plot(timeBinCenters,meanPhaseWidthsErrPerTimeBin,'ko')
hold on
shadedErrorBar(timeBinCenters,meanPhaseWidthsErrPerTimeBin,phaseErrSEM)
title('Linear variablilty model error vs time')
xlabel('Time in field (frac)')
ylabel('Time distribution actual - linear predicted width')
ylim([-60 60])
hold on
plot(xlim,[0 0],'k--','LineWidth',5)

xlim([0 1])
setFigFontTo(18)
maxFigHalfWidth

saveas(gcf,'popAvgTimeVariabilityVsPhase.png')

%%

data.allFieldTimeWidthSlopes

setTightSubplots_Spacious
close all
figure; subplot(2,1,1); 
%offsetDiffsEdges=linspace(-15,15,50);
rEdges=linspace(-1,1,41);
slopeEdges=linspace(-360,360,41);

rBins=edgesToBins(rEdges);
slopeBins=edgesToBins(slopeEdges);

Nr=histcounts(data.allFieldTimeWidthRs,rEdges); 

plot(rBins,smooth(Nr),'r','LineWidth',4)
hold on;
ylim([0 max(smooth(Nr))])
plot([0 0],ylim,'k--','LineWidth',5)
box off
%xlabel('High speed to low speed offset difference (deg)')
xlabel('Width change vs time R')
ylabel('Field count')
xlim([-1 1])
subplot(2,1,2); 
Ns=histcounts(data.allFieldTimeWidthSlopes,slopeEdges); 

plot(slopeBins,smooth(Ns),'r','LineWidth',4)
ylim([0 max(smooth(Ns))])
hold on;
plot([0 0],ylim,'k--','LineWidth',5)
box off
xlabel('Width change per time (deg/field)')
ylabel('Field count')
xlim([-360 360])

setFigFontTo(18)


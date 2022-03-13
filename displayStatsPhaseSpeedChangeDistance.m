data=load('phaseDistanceAndSpeedChangeStatsStrongPrecessing.mat');
statsTable=data.statsTable;
close all

%setTightSubplots

rhoSpeedChange=statsTable{:,1};
slopeSpeedChange=statsTable{:,2};
pSpeedChange=statsTable{:,3};

rhoDist=statsTable{:,4};
slopeDist=statsTable{:,5};
pDist=statsTable{:,6};


rhoEdges=-1:0.1:1;
rhoBins=edgesToBins(rhoEdges);

pEdges=0:0.1:1;
pBins=edgesToBins(pEdges);


alphaLevel=0.05;
sigCorrsSpeedChange=pSpeedChange<alphaLevel;
sigCorrsDist=pDist<alphaLevel;
%rhoSpeedChange=rhoSpeedChange(sigCorrsSpeedChange);
%rhoDist=rhoDist(sigCorrsDist);

%slopeSpeedChange=slopeSpeedChange(sigCorrsSpeedChange);
%slopeDist=slopeDist(sigCorrsDist);



[rhoSpeedChangeN,~]=histcounts((rhoSpeedChange),rhoEdges);
[rhoDistN,~]=histcounts((rhoDist),rhoEdges);

[pSpeedChangeN,~]=histcounts((pSpeedChange),pEdges);
[pDistN,~]=histcounts((pDist),pEdges);

%[pSpeedChangeN,~]=histcounts((pSpeedChange),pEdges);
%[pDistN,~]=histcounts((pDist),pEdges);


%rhoSpeedChangeDensity=rhoSpeedChangeN/sum(rhoSpeedChangeN);
%rhoDistDensity=rhoDistN/sum(rhoDistN);
figure

subplot(1,2,2)
p1=plot(rhoBins,rhoSpeedChangeN/sum(rhoSpeedChangeN),'r','LineWidth',3)
hold on
p2=plot(rhoBins,rhoDistN/sum(rhoDistN),'k','LineWidth',3)

%plot([meanRhoSpeedChange meanRhoSpeedChange],ylim,'r--','LineWidth',5)
%plot([meanRhoDist meanRhoDist],ylim,'k--','LineWidth',5)
xlabel('Circular-linear correlation coefficient')
ylabel('Fraction of place cells')
ylim([0 max(rhoDistN/sum(rhoDistN))])
plot([0 0],ylim,'--','LineWidth',5,'color',getGrayRGB)
title('Strongly precessing place cells correlation \rho')

legend([p1 p2],{'Phase vs change in speed','Phase vs distance'})
subplot(1,2,1)

p1=plot(pBins,pSpeedChangeN/sum(pSpeedChangeN),'r','LineWidth',3)
hold on
p2=plot(pBins,pDistN/sum(pDistN),'k','LineWidth',3)

%plot([meanRhoSpeedChange meanRhoSpeedChange],ylim,'r--','LineWidth',5)
%plot([meanRhoDist meanRhoDist],ylim,'k--','LineWidth',5)
xlabel('Circular-linear correlation p-value')
ylabel('Fraction of place cells')

legend([p1 p2],{'Phase vs change in speed','Phase vs distance'})


title('Strongly precessing place cells correlation p-value')

setFigFontTo(18)
maxFig
saveas(gcf,'PhaseDistanceSpeedChangeCorrelationStats.tif')
%%
figure
subplot(2,1,1)
maxSlopeVal=10000;
slopeEdges=-1000:100:1000;

slopeSpeedChange(abs(slopeSpeedChange)>maxSlopeVal)=[];
slopeDist(abs(slopeDist)>maxSlopeVal)=[];

histogram(slopeSpeedChange,slopeEdges)
subplot(2,1,2)
slopeEdges=-200:10:200;
histogram(slopeDist,slopeEdges)

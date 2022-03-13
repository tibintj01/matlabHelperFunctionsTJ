%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load unit struct file names
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;  clc

%NOT RELOADING LARGE FILE, TAKE OUT TO RELOAD
%clear all;
clearvars -except data spikeDataPerField

%%

%load('allDistTimePhaseTriples1227Fields.mat')
%load('allDistTimePhaseTriples1227FieldsWithMinSpeed.mat')

load('allDistTimePhaseTriples1227FieldsWithMaxAcc.mat')
close all
input.data3cols=allDistTimePhaseTriples;


minNumPts=10;
minNumPts=30;
minNumPts=50;
minNumPts=50;

numBins=50;
numBins=35;
%numBins=20;
useLogBins=1;
useLogBins=0;
%useLogBins=0;

%useLogTimeBins=1;
useLogTimeBins=0;

input.xRad=0.01;
input.yRad=0.01;

input.xRad=0.025;
input.yRad=0.025;

   maxNegExp=1; %base 10
   maxNegExp=7; %base 1.4
   
   timeBase=1.4;
   distBase=1.4;
      %distBase=1.3;
      
if(useLogBins)


    %input.xEdges=10.^(-linspace(0,maxNegExp,numBins+1));
    %input.yEdges=10.^(-linspace(0,maxNegExp,numBins+1));
    
    input.xEdges=timeBase.^(-linspace(0,maxNegExp,numBins+1));
    input.yEdges=distBase.^(-linspace(0,maxNegExp,numBins+1));
else

    %input.xEdges=linspace(0.05,1,numBins+1);
    %input.yEdges=linspace(0.05,1,numBins+1);
    
    input.xEdges=linspace(0.05,0.95,numBins+1);
    input.yEdges=linspace(0.05,0.95,numBins+1);
end

if(useLogTimeBins)
   
      input.yEdges=timeBase.^(-linspace(0,maxNegExp,numBins+1));
      input.xEdges=linspace(1,timeBase^-maxNegExp,numBins+1);
end
input.desiredStat='circMean';
%input.desiredStat='circPeakPhase';

[xCenters,yCenters,allDistTimePhaseMeansSmooth,dataCountPerBin] =makeHeatMapOf3rdVarUniv(input);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot heatmap
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;

allDistTimePhaseMeansSmoothDisp=allDistTimePhaseMeansSmooth;

allDistTimePhaseMeansSmoothDisp(dataCountPerBin<minNumPts)=NaN;
omarPcolor(xCenters,yCenters,allDistTimePhaseMeansSmoothDisp)
colormap(jet)
cb=colorbar

xlabel('Distance in field (frac)')
ylabel('Time in field (frac)')
ylabel(cb,'Circ mean phase (\circ)')

if(useLogBins)
    title('log-log space-time phase precession surface')
    xlabel(sprintf('log_{%.2f}(distance)',distBase))
    ylabel(sprintf('log_{%.2f}(time)',timeBase))
  
end

title(sprintf('All dist-time-phase triplets, HC3 (n=%d fields)',usedFieldCount))

axis square
maxFig
setFigFontTo(16)


%{
if(useLogBins)
    set(gca,'xscale','log')
    set(gca,'yscale','log')
end
if(useLogTimeBins)
    set(gca,'yscale','log')
end
%}


if(useLogTimeBins)
    saveas(gcf,'hc3SpaceLogTimePhase.png')
elseif(useLogBins)
    saveas(gcf,'hc3LogSpaceLogTimePhase.png')
else
    saveas(gcf,'hc3SpaceTimePhase.png')
end

figure;
[distCenters, timeCenters, phaseValues] = prepareSurfaceData(xCenters, yCenters, allDistTimePhaseMeansSmoothDisp);

[fitobject,goodnessOfFit,output] = fit([distCenters,timeCenters],phaseValues,'poly11','Robust','LAR') ;
figure;
plot(fitobject,[distCenters,timeCenters],phaseValues)
colormap(jet)
cb=colorbar
axis square
ylabel(cb,'Circ mean phase (\circ)')

%[b,stats] = robustfit([distCenters(:) timeCenters(:)] ,phaseValues)

title('space-time phase precession surface')
    xlabel('distance')
    ylabel('time')
    zlabel('Theta phase (\circ)')
    
     title({'space-time phase precession surface',...
        sprintf('Best fit plane: theta phase = %.2f%s + (%.2f%s)*(dist) + (%.2f%s)*(time)%s',fitobject.p00,char(176),fitobject.p10,char(176), fitobject.p01,char(176))...
        sprintf('mean absolute error: %.2f%s',mean(abs(output.residuals)),char(176))})
    
if(useLogBins)
    title({'log-log space-time phase precession surface',...
        sprintf('Best fit plane: theta phase = %.2f%s + (%.2f%s)*log_{%.1f}(dist) + (%.2f%s)*log_{%.1f}(time)',fitobject.p00,char(176),fitobject.p10,char(176),distBase, fitobject.p01,char(176),timeBase)...
        sprintf('mean absolute error: %.2f%s',mean(abs(output.residuals)),char(176))})
    xlabel(sprintf('log_{%.2f}(distance)',distBase))
    ylabel(sprintf('log_{%.2f}(time)',timeBase))
    zlabel('Theta phase (\circ)')
end

if(useLogTimeBins)
     title('distance-log(time) phase precession surface')
    xlabel('distance in field')
    ylabel(sprintf('log_{%.2f}(time)',timeBase))
    zlabel('Theta phase (\circ)')
end

xlim([min(distCenters) max(distCenters)])
ylim([min(timeCenters) max(timeCenters)])
zlim([50 360])
setFigFontTo(16)
maxFig

%omarPcolor(xCenters,yCenters,dataCountPerBin)
%colormap(parula)
%colorbar
%%
%%%%%%%%
%plot cross sections
%xCenters,yCenters,allDistTimePhaseMeansSmooth
%close all
distBinWidth=0.1;
timeBinWidth=0.1;

%distBinWidth=1/15;
%timeBinWidth=1/15;

%distBinWidth=1/20;
%timeBinWidth=1/20;

%distBinWidth=1/30;
%timeBinWidth=1/30;


%distBinWidth=1/50;
%timeBinWidth=1/50;



xLabel='Distance in field (norm)';
yLabel='Time in field (norm)';
zLabel='Avg theta phase (deg)';
plotCrossSections(xCenters,yCenters,allDistTimePhaseMeansSmoothDisp,...
    distBinWidth,timeBinWidth,xLabel,yLabel,zLabel)
subplot(1,2,1)
title(sprintf('%s vs %s, holding %s constant','Phase','time','distance'))

xlim([min(input.yEdges) max(input.yEdges)])
ylim([80 340])

subplot(1,2,2)
title(sprintf('%s vs %s, holding %s constant','Phase','distance','time'))
xlim([min(input.xEdges) max(input.xEdges)])
ylim([80 340])

uberTitle(sprintf('n=%d place fields',usedFieldCount))
setFigFontTo(18)
    maxFig
    
    saveas(gcf,'phaseVsTimeAndDistWithControl.png')


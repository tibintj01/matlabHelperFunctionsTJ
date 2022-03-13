close all; clear all; clc
dataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles';

filePaths=getFilePathsRegex(dataDir,'*mat');

allPhases=[];
allSpeeds=[];
for i=1:length(filePaths)
    i
    currFilePath=filePaths{i};
    data=load(currFilePath);
    try
        data.rightSpikesRefPhases;
        data.leftSpikesRefPhases;
    catch
        continue
    end
    
    allSpikePhases=[data.rightSpikesRefPhases(:); data.leftSpikesRefPhases(:)];
    
    allSpikeTimes=[data.rightSpikeTimes(:); data.leftSpikeTimes(:)];
    
    allSpikeSpeeds=abs(interp1(data.positionTimeAxis,data.filledSignedSpeedPerTime,allSpikeTimes));
    
    allPhases=[allPhases;allSpikePhases(:)];
    allSpeeds=[allSpeeds;allSpikeSpeeds(:)]; 

end
%%
minLowSpeed=0.05;
maxLowSpeed=0.1;

minHighSpeed=0.6;
minHighSpeed=0.8;

maxHighSpeed=1.2;

lowSpeedIdxes=allSpeeds>minLowSpeed & allSpeeds<maxLowSpeed;
highSpeedIdxes=allSpeeds>minHighSpeed & allSpeeds<maxHighSpeed;

phaseEdges=linspace(0,360,300);

smoothWind=10;
smoothWind=15;
%smoothWind=25;
smoothWind=30;

%yyaxis left
Nlow=histcounts(allPhases(lowSpeedIdxes),phaseEdges);
%{
Pl=smooth(Nl/sum(Nl),smoothWind);
phaseBins=edgesToBins(phaseEdges);
plot(phaseBins,Pl,'LineWidth',3)
%}

%yyaxis right
Nhigh=histcounts(allPhases(highSpeedIdxes),phaseEdges);
%{
Ph=smooth(Nh/sum(Nh),smoothWind);
plot(phaseBins,Ph,'LineWidth',3)
xlim([0 360])
%}


pPhaseSlow=(Nlow/sum(Nlow(:)));
pPhaseFast=(Nhigh/sum(Nhigh(:)));

pPhaseSlow=[pPhaseSlow(:); pPhaseSlow(:); pPhaseSlow(:); pPhaseSlow(:)];
pPhaseFast=[pPhaseFast(:); pPhaseFast(:); pPhaseFast(:); pPhaseFast(:)];



pPhaseSlow=smooth(pPhaseSlow,smoothWind);
pPhaseFast=smooth(pPhaseFast,smoothWind);


degreesAxis=[edgesToBins(phaseEdges)-360, edgesToBins(phaseEdges), edgesToBins(phaseEdges)+360 , edgesToBins(phaseEdges)+360*2];
figure;
pSlow=plot(degreesAxis,pPhaseSlow,'b-','LineWidth',5)
hold on

pFast=plot(degreesAxis,pPhaseFast,'r-','LineWidth',5)
xlim([-360 720])

plot([0 0],ylim,'k--','LineWidth',5)
plot([0 0]-360,ylim,'k--','LineWidth',5)
plot([0 0]+360,ylim,'k--','LineWidth',5)
plot([0 0]+360*2,ylim,'k--','LineWidth',5)

xlabel('Theta phase of place cell spiking')
ylabel('Probability')
title('All place fields in population')
legend([pSlow pFast],sprintf('running speed < %.2f m/s',maxLowSpeed),sprintf('running speed > %.2f m/s',minHighSpeed))

%xlim([0 360])
%xlim([60 300])
setFigFontTo(32)
maxFig

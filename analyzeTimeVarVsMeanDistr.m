close all; clear all; clc
processedDataDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData';

saveDataPath=fullfile(processedDataDir,'phaseDistrWidthVsTimeInFieldData.mat');

data=load(saveDataPath);

allFieldTimeDistrMeans=data.allFieldTimeDistrMeans;
allFieldTimeDistrWidths=data.allFieldTimeDistrWidths;

phaseRepInvariantTimeRep=1;



if(phaseRepInvariantTimeRep)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %phase binned, mean labeled time values
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    numTimeRepresentBins=7;
     numTimeRepresentBins=12
       %numTimeRepresentBins=11
      %numTimeRepresentBins=12


    timeRepresentEdges=linspace(0,360,numTimeRepresentBins+1);
    phaseBinCenters=edgesToBins(timeRepresentEdges);
    phaseBinWidth=median(diff(phaseBinCenters));
    %allTimes=data
    allUsedFieldsPhases=data.allUsedFieldsPhases;
    allUsedFieldsTimes=data.allUsedFieldsTimes;

    figure;
    timeRepColors=jet(numTimeRepresentBins);
    numActualTimeBins=50;
      numActualTimeBins=40;
    %numActualTimeBins=100;
    cutOffTime=0;
    actualTimeEdges=linspace(0,1,numActualTimeBins+1);
     %actualTimeEdges=linspace(-0.2,1,numActualTimeBins+1);
    actualTimeBins=edgesToBins(actualTimeEdges);

    for ti=1:numTimeRepresentBins
       currStartPhase=phaseBinCenters(ti)-phaseBinWidth/2;
       currEndPhase=phaseBinCenters(ti)+phaseBinWidth/2;
       
       currPhaseBinIdxes=allUsedFieldsPhases>currStartPhase & allUsedFieldsPhases<=currEndPhase;
       
       currPhaseBinTimes=allUsedFieldsTimes(currPhaseBinIdxes);
       

       Nt=histcounts(currPhaseBinTimes,actualTimeEdges);
       smoothIdx=5;
       currTimeRepDistr=smooth(Nt,smoothIdx)/sum(smooth(Nt,smoothIdx));
       
       currTimeRepDistr=currTimeRepDistr/max(currTimeRepDistr(:));


       if(ti>2 & ti<numTimeRepresentBins-1 )
         plot(actualTimeBins,currTimeRepDistr,'LineWidth',4,'Color',timeRepColors(ti,:))
       end
       hold on

    end
    xlim([0 0.9])
    title('peak normalized')
    disp('')
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %phase binned, mean labeled time values
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    numTimeRepresentBins=7;
    numTimeRepresentBins=10;
    numTimeRepresentBins=20;

    timeRepresentEdges=linspace(-0.1,1,numTimeRepresentBins);
    phaseBinnedMeanLabels=data.allFieldTimeMeansByPhaseBins;
    phaseBinnedTimeValues=data.allFieldTimeValsByPhaseBins;

    timeRepBinLabels=discretize(phaseBinnedMeanLabels,timeRepresentEdges);

    figure;
    timeRepColors=jet(numTimeRepresentBins);
    numActualTimeBins=50;
    %numActualTimeBins=100;
    cutOffTime=0;
    actualTimeEdges=linspace(-0.1+cutOffTime,1-cutOffTime,numActualTimeBins);
    actualTimeBins=edgesToBins(actualTimeEdges);
    edgeCutOffBins=5;
    for ti=1:numTimeRepresentBins
       currTimeRepBinIdxes=(timeRepBinLabels==ti);

       currTimeRepBinActualTimes=phaseBinnedTimeValues(currTimeRepBinIdxes);

       Nt=histcounts(currTimeRepBinActualTimes,actualTimeEdges);
       smoothIdx=5;
       currTimeRepDistr=smooth(Nt,smoothIdx)/sum(smooth(Nt,smoothIdx));

       if(ti>edgeCutOffBins & ti <numTimeRepresentBins-edgeCutOffBins)
      plot(actualTimeBins,currTimeRepDistr,'LineWidth',4,'Color',timeRepColors(ti,:))
       hold on
       end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tolLev=0.01;

numFieldsUsed=length(data.allFieldTimeWidthRs);

%allFieldLinearModelErrs=data.allFieldLinearModelErrs;

numTimeBins=50;
numTimeBins=75;
numTimeBins=60;
numTimeBins=50;

%meanTimeBinEdges=linspace(0,0.5,numTimeBins+1);
minVarTime=0
%minVarTime=0.1;
maxVarTime=0.4;
meanTimeBinEdges=linspace(minVarTime,maxVarTime,numTimeBins+1);

varTimeBinEdges=linspace(0,1,numTimeBins+1);

fH=figure
subplot(2,2,[1 2])
%plot(allFieldTimeDistrMeans,allFieldTimeDistrWidths,'k.')
getJointDistr(allFieldTimeDistrMeans,allFieldTimeDistrWidths,varTimeBinEdges,meanTimeBinEdges,fH)
xlim([0 1])
ylim([0.1 0.35])
ylim([minVarTime maxVarTime])
%daspect([1 1 1])
axis square
xlabel('Mean of time distribution')
ylabel('Std. dev. of time distribution')

%ylabel(cb,'Joint p

subplot(2,3,4)
rEdges=linspace(-1,1,101);
rBins=edgesToBins(rEdges);

slopeEdges=linspace(-1,1,101);
slopeBins=edgesToBins(slopeEdges);

offsetEdges=linspace(-1,1,501);
offsetBins=edgesToBins(offsetEdges);

Ns=histcounts(data.allFieldTimeWidthSlopes,slopeEdges);

plot(slopeBins,smooth(Ns/sum(Ns(:)),3),'r-','LineWidth',4)
axis square
hold on

plot([0 0],ylim,'k--','LineWidth',4)
ylim([-Inf Inf])
xlabel('mean vs std of time distributions, slope')

ylabel('Probability')

subplot(2,3,5)
Nr=histcounts(data.allFieldTimeWidthRs,rEdges);
plot(rBins,smooth(Nr/sum(Nr(:)),3),'r-','LineWidth',4)
axis square
hold on

plot([0 0],ylim,'k--','LineWidth',4)
ylim([-Inf Inf])

xlabel('mean vs std of time distributions, correlation coefficient')
ylabel('Probability')


subplot(2,3,6)
No=histcounts(data.allFieldTimeWidthOffsets,offsetEdges);
plot(offsetBins,smooth(No/sum(No(:)),3),'r-','LineWidth',4)
axis square
hold on

plot([0 0],ylim,'k--','LineWidth',4)
ylim([-Inf Inf])

xlabel('mean vs std of time distributions, offset')
ylabel('Probability')

maxFig
setFigFontTo(24)

uberTitle({'Scale invariant representation of time (linear relationship between mean and std)', sprintf('within theta cycles (n=%d place fields)',numFieldsUsed)},18)

%meanPhaseWidthsPerTimeBin=sqrt(meanPhaseWidthsPerTimeBin);
setTightSubplots


saveas(gcf,'scaleInvariantTimeCoding.png')




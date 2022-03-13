close all; clear all; clc
processedDataDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData';

saveDataPath=fullfile(processedDataDir,'phaseDistrWidthVsTimeInFieldData.mat');

data=load(saveDataPath);

allFieldTimeBinCenters=data.allFieldTimeBinCenters;
timeBinCenters=unique(allFieldTimeBinCenters);

allFieldDistrWidths=data.allFieldDistrWidths;

allFieldDistrWidths(allFieldDistrWidths>360)=NaN;

tolLev=0.01;

meanPhaseWidthsPerTimeBin=NaN(length(timeBinCenters),1);
meanPhaseWidthsErrPerTimeBin=NaN(length(timeBinCenters),1);
meanPhaseWidthsQuadErrPerTimeBin=NaN(length(timeBinCenters),1);

phaseErrSEM=NaN(length(timeBinCenters),1);
phaseQuadErrSEM=NaN(length(timeBinCenters),1);
phaseSEM=NaN(length(timeBinCenters),1);

allFieldLinearModelErrs=data.allFieldLinearModelErrs;
allFieldQuadModelErrs=data.allFieldQuadModelErrs;

figure
for i=1:length(timeBinCenters)
    

    currTimeBinIdxes=(abs(allFieldTimeBinCenters-timeBinCenters(i))<tolLev);
    
    currTimeBinLinModelErrs=allFieldLinearModelErrs(currTimeBinIdxes);
    %currTimeBinQuadModelErrs=allFieldQuadModelErrs(currTimeBinIdxes);
    
    currTimeBinWidths=allFieldDistrWidths(currTimeBinIdxes);
    %currTimeBinWidths=sqrt(currTimeBinWidths-120);
    %meanPhaseWidthsPerTimeBin(i)=circMeanDeg(currTimeBinWidths);
    %phaseSEM(i)=1/circ_kappa(currTimeBinWidths);
    meanPhaseWidthsPerTimeBin(i)=nanmean(currTimeBinWidths);
    meanPhaseWidthsErrPerTimeBin(i)=nanmean(currTimeBinLinModelErrs);
    meanPhaseWidthsQuadErrPerTimeBin(i)=nanmean(currTimeBinQuadModelErrs);
    
    phaseErrSEM(i)=getSEMacrossRows(currTimeBinLinModelErrs(:));
    phaseQuadErrSEM(i)=getSEMacrossRows(currTimeBinQuadModelErrs(:));
    
    phaseSEM(i)=getSEMacrossRows(currTimeBinWidths(:));
    
        subplot(4,4,i)
        histogram(currTimeBinWidths,linspace(90,180,30))
    
end

%meanPhaseWidthsPerTimeBin=sqrt(meanPhaseWidthsPerTimeBin);
setTightSubplots_Spacious

figure;
subplot(3,1,1)
plot(timeBinCenters,meanPhaseWidthsPerTimeBin,'ko')
hold on
shadedErrorBar(timeBinCenters,meanPhaseWidthsPerTimeBin,phaseSEM)
title('Phase distribution width vs time in field')

xlabel('Time in field (frac)')
ylabel('Phase distribution width (deg)')

subplot(3,1,2)
plot(timeBinCenters,meanPhaseWidthsErrPerTimeBin,'ko')
shadedErrorBar(timeBinCenters,meanPhaseWidthsErrPerTimeBin,phaseErrSEM)
title('Linear variablilty model error vs time in field')
xlabel('Time in field (frac)')
ylabel('Phase distribution actual - linear predicted width (deg)')
ylim([-20 20])

subplot(3,1,3)
plot(timeBinCenters,meanPhaseWidthsQuadErrPerTimeBin,'ko')
shadedErrorBar(timeBinCenters,meanPhaseWidthsQuadErrPerTimeBin,phaseQuadErrSEM)
title('Quadratic variablilty model error vs time in field')
xlabel('Time in field (frac)')
ylabel('Phase distribution actual - quadratic predicted width (deg)')
ylim([-20 20])

setFigFontTo(12)


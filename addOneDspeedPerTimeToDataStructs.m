close all; clear all; clc
dataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles';
showPlots=0;
debugTimeFilePath='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles/Achilles_11012013_6-12Theta_Unit1024Info.mat';
debugTimeFilePath='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles/Achilles_11012013_6-12Theta_Unit1006Info.mat';
debugTimeFilePath='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles/Achilles_11012013_6-12Theta_Unit1024Info.mat';
debugTimeFilePath='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles/Achilles_11012013_6-12Theta_Unit1038Info.mat';
debugTimeFilePath='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles/Achilles_11012013_6-12Theta_Unit1114Info.mat';
debugTimeFilePath='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles/Cicero_09102014_6-12Theta_Unit929Info.mat';
debugTimeFilePath='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles/Achilles_11012013_6-12Theta_Unit1006Info.mat';
debugTimeFilePath='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles/Achilles_11012013_6-12Theta_Unit102Info.mat';
debugTimeFilePath='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles/Achilles_10252013_6-12Theta_Unit1115Info.mat';

filePaths=getFilePathsRegex(dataDir,'*mat');

for i=1:length(filePaths)

    currFilePath=filePaths{i};
    i/length(filePaths)

    data=load(currFilePath);

    oneDLoc=data.unitInfo.positionPerTime.OneDLocation(:);
    oneDtimeAxis=data.unitInfo.positionPerTime.TimeStamps(:);

    approxPosTimeStep=median(diff(oneDtimeAxis));

    polynomialFitOrder=1; %fit lines over acc time width
    speedEstTimeWidth=0.125; %sec %~5 point linear fit estimation of speed
      %speedEstTimeWidth=0.075; %sec %~5 point linear fit estimation of speed
    supportWindLength=ceil(speedEstTimeWidth/approxPosTimeStep);
    oneDspeedPerTimeStep=movingslope(oneDLoc,supportWindLength,polynomialFitOrder,approxPosTimeStep);
       
    absOneDspeedPerTimeStep=abs(oneDspeedPerTimeStep);
 
    %{
    figure
    plot(oneDtimeAxis,absOneDspeedPerTimeStep,'k-')
    hold on
    plot(oneDtimeAxis,oneDLoc,'b-')
%}

    save(currFilePath,'oneDtimeAxis','oneDLoc','oneDspeedPerTimeStep','absOneDspeedPerTimeStep','-append')

    
end
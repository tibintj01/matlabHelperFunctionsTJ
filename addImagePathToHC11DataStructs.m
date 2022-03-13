close all; clear all; clc
dataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles';

imageDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/placeCellQualityImages/singleDirectionGoodPrecessionImages';
filePaths=getFilePathsRegex(dataDir,'*mat');

%for fi=1:length(filePaths)
for fi=1:1
       
    currFilePath=filePaths{fi};
    
    currFileName=getFileNameFromPath(currFilePath);
    fileBaseName=currFileName(1:(end-4));
    
    imagePaths=getFilePathsRegex(imageDir,sprintf('%s*.png',fileBaseName));
    
    save(currFilePath,'imagePaths','-append')
end
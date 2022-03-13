clear all; close all; clc
imageDir='/Users/tibinjohn/thetaSeq/code/spaceTimeProcessedDir/manuallySortedPlaceCellsAbove10Hz/strong/tifs';
dataDir='/Users/tibinjohn/thetaSeq/code/spaceTimeProcessedDir/manuallySortedPlaceCellsAbove10Hz/strong';
saveDir='/Users/tibinjohn/thetaSeq/code/PhaseAndDistVsPhaseAndSpeedChange';

dataFilePaths=getFilePathsRegex(dataDir,'*mat');
imageFilePaths=getFilePathsRegex(imageDir,'*tif');

touchDir(saveDir)

statsTable=table(NaN,NaN,NaN,NaN,NaN,NaN,'VariableNames',{'circLinRhoSpeed','circLinSlopeSpeed','circLinPSpeed', 'circLinRhoDist','circLinSlopeDist','circLinPDist'})
for fi=1:length(dataFilePaths)
    currFilePath=dataFilePaths{fi};
    currImgPath=imageFilePaths{fi};
    currData=load(currFilePath)
    
    currStartCm=currData.manualFieldStartCm;
    currEndCm=currData.manualFieldEndCm;
    manualPrecessionStartPhaseDeg=currData.manualPrecessionStartPhaseDeg;
    
    if(currEndCm<=115)
        %continue
    end


    fileSaveRootName=getFileNameFromPath(currFilePath);
    saveRootName=fileSaveRootName(1:(end-4));

    %getSpeedChangeVsPosition(currData)
    [circLinRhoSpeed,circLinSlopeSpeed,circLinPSpeed, circLinRhoDist,circLinSlopeDist,circLinPDist]=getSpeedChangeVsPositionSaturatingFit(currData);

        currTableRow=table(circLinRhoSpeed,circLinSlopeSpeed,circLinPSpeed, circLinRhoDist,circLinSlopeDist,circLinPDist);
  
     statsTable=[statsTable;   currTableRow];
    %saveas(gcf,fullfile(saveDir, [saveRootName '.tif']))
    

  close all
end
statsTable(1,:)=[];
save('phaseDistanceAndSpeedChangeStatsStrongPrecessing.mat','statsTable')

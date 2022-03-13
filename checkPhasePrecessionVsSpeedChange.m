imageDir='/Users/tibinjohn/thetaSeq/code/spaceTimeProcessedDir/manuallySortedPlaceCellsAbove10Hz/strong/tifs';
dataDir='/Users/tibinjohn/thetaSeq/code/spaceTimeProcessedDir/manuallySortedPlaceCellsAbove10Hz/strong';

dataFilePaths=getFilePathsRegex(dataDir,'*mat');
imageFilePaths=getFilePathsRegex(imageDir,'*tif');

centerFieldsDir=fullfile(imageDir,'centerFields');
earlyFieldDir=fullfile(imageDir,'earlyFields');
lateFieldDir=fullfile(imageDir,'lateFields');

touchDir(centerFieldsDir)
touchDir(earlyFieldDir)
touchDir(lateFieldDir)

close all
middleStartCm=35;
middleEndCm=85;
middleCm=60;
endCm=100;

for fi=1:length(dataFilePaths)
    currFilePath=dataFilePaths{fi};
    currImgPath=imageFilePaths{fi};
    currData=load(currFilePath)
    
    currStartCm=currData.manualFieldStartCm;
    currEndCm=currData.manualFieldEndCm;
    manualPrecessionStartPhaseDeg=currData.manualPrecessionStartPhaseDeg;


    fileSaveRootName=getFileNameFromPath(currFilePath);
    saveRootName=fileSaveRootName(1:(end-4));
    if(currData.spaceTimePhaseInfo.dirFlag<0)
       di=1;
    else
       di=0;
    end

    rateCOM=getCOM(currData.spaceTimePhaseInfo.placeField);

    plotPrecessionAndSpeed

    
    if(currStartCm>=middleStartCm && currEndCm < middleEndCm)
        %figure
        %imshow(currImgPath)
	saveas(gcf,[saveRootName 'centerField.tif'])
	 copyfile(currImgPath,centerFieldsDir)
    elseif(currEndCm<middleCm && rateCOM < middleCm)
        copyfile(currImgPath,earlyFieldDir)
	saveas(gcf,[saveRootName 'earlyField.tif'])
    elseif(currEndCm>=endCm && rateCOM > middleCm+10)
        copyfile(currImgPath,lateFieldDir)
	saveas(gcf,[saveRootName 'lateField.tif'])
    end
    disp('here')
  
end

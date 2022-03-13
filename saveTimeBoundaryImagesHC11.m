imageSaveDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/placeCellQualityImages/singleDirectionGoodPrecessionImagesTime';
dataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles';

touchDir(imageSaveDir)

filePaths=getFilePathsRegex(dataDir,'*Info.mat');

for fi=1:length(filePaths)
    currFilePath=filePaths{fi};
    fi
    if(fi==3)
        currUnitInfo=load(currFilePath)
        %continue
    end
    currUnitInfo=load(currFilePath);
    plotTimeFieldAndThetaPhases(currUnitInfo,imageSaveDir);
end


destDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles';
sourceDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/'




  destDataFilePaths=getFilePathsRegex(destDir,'*Info.mat');

for fi=1:length(destDataFilePaths)
	currFilePath=destDataFilePaths{fi};

	currFileName=getFileNameFromPath(currFilePath);
    
    if(~strContainsCircSessionName(currFilePath))
        continue
    end
    updatedFilePath=fullfile(sourceDir,currFileName);
    
    copyfile(updatedFilePath,currFilePath)
    
    
end


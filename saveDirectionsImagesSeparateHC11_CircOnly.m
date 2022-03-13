imageSaveDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/placeCellQualityImages/singleDirectionGoodPrecessionImages/redoCirc';
dataBaseDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/';

touchDir(imageSaveDir)

dirNames={'strongLeftward','strongRightward','strongBothWays'};


for dirIdx=1:length(dirNames)
	currDirName=dirNames{dirIdx};

	currDirFilePaths=getFilePathsRegex(fullfile(dataBaseDir),'*Info.mat');

	for fi=1:length(currDirFilePaths)
	    currFilePath=currDirFilePaths{fi};
	    currUnitInfo=load(currFilePath)
        plotPlaceFieldAndThetaPhases(currUnitInfo,'betterLeftward',imageSaveDir);
        
        %{
		if(strcmp(currDirName,'strongLeftward'))
		    plotPlaceFieldAndThetaPhases(currUnitInfo,'betterLeftward',imageSaveDir);
		elseif(strcmp(currDirName,'strongRightward'))
		    plotPlaceFieldAndThetaPhases(currUnitInfo,'betterRightward',imageSaveDir);
		elseif(strcmp(currDirName,'strongBothWays'))
		     plotPlaceFieldAndThetaPhases(currUnitInfo,'betterLeftward',imageSaveDir);
		     plotPlaceFieldAndThetaPhases(currUnitInfo,'betterRightward',imageSaveDir);
		end
        %}
	end
end

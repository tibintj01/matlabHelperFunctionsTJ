close all; clear all; clc
dataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles';
dataBaseDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz';
dataBaseDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData';

%dirNames={'strongLeftward','strongRightward','strongBothWays'};

%for dirIdx=1:length(dirNames)
%	currDirName=dirNames{dirIdx};

%filePaths=getFilePathsRegex(dataDir,'*mat');
	%filePaths=getFilePathsRegex(dataDir,'*mat');
	%filePaths=getFilePathsRegex(fullfile(dataBaseDir,currDirName),'*Info.mat');
	filePaths=getFilePathsRegex(fullfile(dataBaseDir),'*Info.mat');

	for fi=1:length(filePaths)
	    currFilePath=filePaths{fi};
	    unitData=load(currFilePath);

	    cycleData=load(unitData.unitInfo.cycleInfoFileName);
	    ampPerCycle=cycleData.cycleMaxAmp-cycleData.cycleMinAmp;
	    durPerCycle=cycleData.cycleDurations;
	    minMaxDurPerCycle=cycleData.cycleMinMaxDurations;
	    cycleStartTimes=cycleData.cycleMaxTimes;
	    

	    for di=1:2
		    if(di==1)
			allSpikeTimes=unitData.rightSpikeTimes;
		    elseif(di==2)
			allSpikeTimes=unitData.leftSpikeTimes;
		    end

		    cycleAmpPerSpike=NaN(size(allSpikeTimes));
		    cycleDurPerSpike=NaN(size(allSpikeTimes));
		    cycleMinMaxDurPerSpike=NaN(size(allSpikeTimes));

		    for si=1:length(allSpikeTimes)
		    currSpikeTime=allSpikeTimes(si);
		    [~,closestCycleIdx]=min(abs(currSpikeTime-cycleStartTimes));

		    latestCycleStartTime=cycleStartTimes(closestCycleIdx);
		    if(currSpikeTime<latestCycleStartTime && closestCycleIdx>1)
			latestCycleStartTime=cycleStartTimes(closestCycleIdx-1);
			closestCycleIdx=closestCycleIdx-1;
		    end

		    cycleAmpPerSpike(si)=ampPerCycle(closestCycleIdx);
		    cycleDurPerSpike(si)=durPerCycle(closestCycleIdx);
		    cycleMinMaxDurPerSpike(si)=minMaxDurPerCycle(closestCycleIdx);
		end

		if(di==1)
		    cycleAmpPerRightSpike=cycleAmpPerSpike;
		    cycleDurPerRightSpike=cycleDurPerSpike;
		    cycleMinMaxDurPerRightSpike=cycleMinMaxDurPerSpike;
		elseif(di==2)
		    cycleAmpPerLeftSpike=cycleAmpPerSpike;
		    cycleDurPerLeftSpike=cycleDurPerSpike;
		    cycleMinMaxDurPerLeftSpike=cycleMinMaxDurPerSpike;
		end

	     end
	    save(currFilePath,'cycleAmpPerLeftSpike','cycleDurPerLeftSpike','cycleMinMaxDurPerLeftSpike',...
		'cycleAmpPerRightSpike','cycleDurPerRightSpike','cycleMinMaxDurPerRightSpike','-append')
	    dataCheck=load(currFilePath)
	end

%end

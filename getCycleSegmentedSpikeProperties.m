function []=saveCycleSegmentedSpikeProperties(ptDir,sessionNum,ch,singleCycleType)

	cycleSegmentedSpikeProperties.ch=ch;
	cycleSegmentedSpikeProperties.ptDir=ptDir;
	cycleSegmentedSpikeProperties.sessionNum=sessionNum;
	cycleSegmentedSpikeProperties.singleCycleType=singleCycleType;	

	cellPropDir=sprintf('/nfs/turbo/lsa-ojahmed/processedHumanData/%s/sessionID-%d/cellProperties-MatFiles',ptDir,sessionNum);

	singleCycleMatDirPath=sprintf('/nfs/turbo/lsa-ojahmed/processedHumanData/%s/sessionID-%d/singleCycleProperties-MatFiles/%s',ptDir,sessionNum,singleCycleType);

	if(ch<10)
		chStr=['0' num2str(ch)];
	else
		chStr=num2str(ch);
	end

	currChCellFileNames=getRegexFilePaths(cellPropDir,sprintf('*%s*_cell*',chStr));
	currChSingleCycleFileName=getRegexFilePaths(singleCycleMatDirPath,sprintf('*Ch%s.mat',chStr));

	singleCycleInfo=load(currChSingleCycleFileName{1});
	
	numCycles=length(singleCycleInfo.asymData.fastStartTimes);
	%fds
	numCells=length(currChCellFileNames);

	cycleSpikeCounts=zeros(numCells,numCycles);
	maxNumSpikesInCycle=20;
	spikePhaseMatrix=-ones(numCells,numCycles,maxNumSpikesInCycle);	
	%loop through cell prop files and populate cycleSegmentedSpikeProperties structure
	for cellNum=1:numCells
		cellPropFile=currChCellFileNames{cellNum};
		cellProps=load(cellPropFile);

		for cycleNum=1:numCycles
			cycleStart=singleCycleInfo.asymData.fastStartTimes(cycleNum);
			cycleEnd=singleCycleInfo.asymData.fastEndTimes(cycleNum);
			cycleSpikeCounts(cellNum,cycleNum)=sum(cellProps.spikeTimes>=cycleStart & cellProps.spikeTimes<=cycleEnd);
			if(cycleSpikeCounts(cellNum,cycleNum)==0)
				continue;
			end

			withinCycleSpikeTimes=cellProps.spikeTimes(cellProps.spikeTimes>=cycleStart & cellProps.spikeTimes<=cycleEnd)-cycleStart;

			withinCycleSpikePhases=withinCycleSpikeTimes/(cycleEnd-cycleStart)*360; %0deg is start, 360 is end
			spikePhaseMatrix(cellNum,cycleNum,1:cycleSpikeCounts(cellNum,cycleNum))=withinCycleSpikePhases;

		end
		%fds
	end

	 cycleSegmentedSpikeProperties.cycleSpikeCounts=cycleSpikeCounts;
       	cycleSegmentedSpikeProperties.spikePhaseMatrix=spikePhaseMatrix;
	cycleSegmentedSpikeProperties.currChCellFileNames=currChCellFileNames;

    save(sprintf('singleCycle%sDurationNormed-Ch%sSpikeProps.mat',singleCycleType,chStr),'cycleSegmentedSpikeProperties')
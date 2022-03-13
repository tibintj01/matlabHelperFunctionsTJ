%singleStatsDir='/nfs/turbo/lsa-ojahmed/mg49_singleCycleStats_Par';
singleStatsDir='/nfs/turbo/lsa-ojahmed/mg49_singleCycleBroadBandStats_Par';

freqBandName='Alpha';
singleCycleMatFilePaths=getRegexFilePaths(singleStatsDir,sprintf('%s*mat',freqBandName));

numChannels=length(singleCycleMatFilePaths);
%numChannels=20
numChannels=6
numChannels=1
%numChannels=50

%200 ms window
%halfCycleWindow=0.1;
halfCycleWindow=0.075;
%halfCycleWindow=0.05;

zscoreCycles=1;

if(zscoreCycles)
    zscoreStr='zScoredCycles';
else
    zscoreStr='notZscoredCycles';
end

decFs=30000/15;

halfCycleIdxLength=halfCycleWindow*(decFs);

disp('calculating size of matrix.....')
tic
numCycles=0;
for i=1:numChannels
	data=load(singleCycleMatFilePaths{i});
	%numCycles=numCycles+length(data.asymData.ripMaxTimes);
	numCycles=numCycles+length(data.asymData.cycleMidTimes);
end
toc

cycleMatrix=zeros(round(2*halfCycleWindow*decFs)+1,numCycles);

timeAxis=(0:(1/decFs):2*halfCycleWindow) - halfCycleWindow;

size(cycleMatrix)
totalNumCyclesSaved=0;
for i=1:numChannels
	fileName=getFileNameFromPath(singleCycleMatFilePaths{i});
	data=load(singleCycleMatFilePaths{i});
	ch=str2num(fileName((end-5):(end-4)));
	[decLFP,decFs]=getDecLFPforCh(ch);

	[lowFastFreq,highFastFreq]=getFastFreqBand(freqBandName);
	%decLFP=filterLFP(decLFP,decFs,data.lowFreq,data.highFreq,2,1);
	decLFP=filterLFP(decLFP,decFs,lowFastFreq,highFastFreq,2,1);
	lfpLength=length(decLFP);
	%numCyclesForCh=length(data.singleCycleInfo.ripMaxTimes);	
	numCyclesForCh=length(data.asymData.fastMidTimes);	
	disp(sprintf('Loading ch %d........',i))
	tic	
	for cycleNum=1:numCyclesForCh
		%cycleStartIdx=(data.asymData.ripMaxPoints(cycleNum)-halfCycleIdxLength);
		%cycleEndIdx=(data.asymData.ripMaxPoints(cycleNum)+halfCycleIdxLength);
		cycleStartIdx=(data.asymData.fastMidInd(cycleNum)-halfCycleIdxLength);
		cycleEndIdx=(data.asymData.fastMidInd(cycleNum)+halfCycleIdxLength);
		%cycleMax=data.asymData.ripMaxAmp(cycleNum);
		cycleMax=data.asymData.fastMidAmp(cycleNum);
		if(cycleStartIdx>0 && cycleEndIdx <=lfpLength)
		    if(zscoreCycles)
			cycleMatrix(:,cycleNum)=zscoreLFP(decLFP(cycleStartIdx:cycleEndIdx));
		    else
			cycleMatrix(:,cycleNum)=decLFP(cycleStartIdx:cycleEndIdx)/cycleMax;
		    end
			totalNumCyclesSaved=totalNumCyclesSaved+1;
		end
	end
	toc
end

%cycleMatrix=subtractMeanCol(cycleMatrix);
figure
plot(timeAxis,cycleMatrix(:,1:100))

%fds
cycleMatrix=cycleMatrix(:,1:totalNumCyclesSaved);

%fds
disp('saving cycle matrix.....')
tic
save(sprintf('cycleMatrixBroadbandZscore-MG49-%s-%s-%0.1f-ms-Window.mat',freqBandName,zscoreStr,halfCycleWindow*2000),'cycleMatrix','-v7.3')

toc

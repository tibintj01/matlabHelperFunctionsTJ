%singleStatsDir='/nfs/turbo/lsa-ojahmed/mg49_singleCycleStats_Par';
%singleStatsDir='/nfs/turbo/lsa-ojahmed/mg49_singleCycleBroadBandStats_Par';

singleStatsDir='/nfs/turbo/lsa-ojahmed/processedHumanData/MG49/sessionID-3/singleCycleProperties-MatFiles/broadbandAlpha';

%get all single cycle property file paths correpsonding to rhythm band
freqBandName='Alpha';
singleCycleMatFilePaths=getRegexFilePaths(singleStatsDir,sprintf('%s*mat',freqBandName));

numChannels=length(singleCycleMatFilePaths);
%numChannels=20
numChannels=6
numChannels=1
%numChannels=50


zscoreCycles=1;

if(zscoreCycles)
    zscoreStr='zScoredCycles';
else
    zscoreStr='notZscoredCycles';
end

%samples/sec
decFs=30000/15;

%cycles/sec
[fLow fHigh]=getFastFreqBand(freqBandName);

%samples/cycle
fixedNumPoints=ceil(decFs/fLow);


disp('calculating size of matrix.....')
tic
numCycles=0;
for i=1:numChannels
	data=load(singleCycleMatFilePaths{i});
	%numCycles=numCycles+length(data.asymData.ripMaxTimes);
	numCycles=numCycles+length(data.asymData.cycleMidTimes);
end
toc

%initialize matrix and time axis for population by data
%cycleMatrix=zeros(round(2*halfCycleWindow*decFs)+1,numCycles);
cycleMatrix=zeros(fixedNumPoints,numCycles);

%timeAxis=(0:(1/decFs):2*halfCycleWindow) - halfCycleWindow;
timeAxis=linspace(0,1,fixedNumPoints);

actualDurations=zeros(numCycles,1);

size(cycleMatrix)

%populate matrix with duration-normalized cycle data, looping over channels
totalNumCyclesSaved=0;
for i=1:numChannels
	disp(sprintf('extracting cycles from ch %d........',i))
	%load single cycle info
	data=load(singleCycleMatFilePaths{i});

	%extract raw decimated LFP corresponding to current channel
	fileName=getFileNameFromPath(singleCycleMatFilePaths{i});
	ch=str2num(fileName((end-5):(end-4)));
	[decLFP,decFs]=getDecLFPforCh(ch);

	if(useBroadband)
		%get filtered LFP in desired (broad) frequency band
		[lowFastFreq,highFastFreq]=getFastFreqBand(freqBandName);
		%decLFP=filterLFP(decLFP,decFs,data.lowFreq,data.highFreq,2,1);
		decLFP=filterLFP(decLFP,decFs,lowFastFreq,highFastFreq,2,1);
	else
		%get filtered LFP in desired (broad) frequency band
		[lowFreq,highFreq]=getFreqBand(freqBandName);
		%decLFP=filterLFP(decLFP,decFs,data.lowFreq,data.highFreq,2,1);
		decLFP=filterLFP(decLFP,decFs,lowFreq,highFreq,2,1);

	end

	lfpLength=length(decLFP);
	%numCyclesForCh=length(data.singleCycleInfo.ripMaxTimes);	
	numCyclesForCh=length(data.asymData.fastMidTimes);	

		
	for cycleNum=1:numCyclesForCh
		%cycleStartIdx=(data.asymData.fastMidInd(cycleNum)-halfCycleIdxLength);
		%cycleEndIdx=(data.asymData.fastMidInd(cycleNum)+halfCycleIdxLength);
		cycleStartIdx=data.asymData.fastStartInd(cycleNum);
		cycleEndIdx=data.asymData.fastEndInd(cycleNum);

		cycleDuration=data.asymData.fastStartEndDuration(cycleNum);

		resamplePadLength=floor((cycleEndIdx-cycleStartIdx)/2);
		
		%cycleMax=data.asymData.ripMaxAmp(cycleNum);
		cycleMax=data.asymData.fastMidAmp(cycleNum);
		%if(cycleStartIdx>1 && cycleEndIdx <=(lfpLength-1))
		if(cycleStartIdx>resamplePadLength && cycleEndIdx <=(lfpLength-resamplePadLength))
		    if(zscoreCycles)
			%startToEndCycle=zscoreLFP(decLFP(cycleStartIdx:cycleEndIdx));
			%resample assumes 0's at either end, take 1 data point as padding
			%startToEndCycle=zscoreLFP(decLFP((cycleStartIdx-1):(cycleEndIdx+1)));
			%startToEndCycle=zscoreLFP(decLFP((cycleStartIdx-1):(cycleEndIdx+1)));
			rawStartToEndCyclePadded=decLFP((cycleStartIdx-resamplePadLength):(cycleEndIdx+resamplePadLength));
			%actualLength=length(startToEndCycle);
			%actualLength=length(rawStartToEndCyclePadded);
			actualLength=length(rawStartToEndCyclePadded)-resamplePadLength*2;
			%desiredLength=fixedNumPoints+resamplePadLength*2;
			desiredLength=fixedNumPoints;
			%disp('resampling cycle...')
			%resampledCycle=resample(startToEndCycle,desiredLength,actualLength); 
			resampledCycle=resample(rawStartToEndCyclePadded,desiredLength,actualLength); 
			resampledResamplePadLength=ceil(resamplePadLength*(desiredLength/actualLength));
			%cycleMatrix(:,cycleNum)=resampledCycle((resamplePadLength+1):(end-resamplePadLength));
			%may introduce TIME JITTER of up to 1 step (0.5 ms here) 
			if(length(resampledCycle((resampledResamplePadLength):(end-resampledResamplePadLength)))==desiredLength)
				cycleMatrix(:,cycleNum)=resampledCycle((resampledResamplePadLength):(end-resampledResamplePadLength));
				%disp('not rounding off...')
			else
				cycleMatrix(:,cycleNum)=resampledCycle((resampledResamplePadLength):(end-resampledResamplePadLength-1));
				%disp('rounding off...')
			end
			cycleMatrix(:,cycleNum)=zscoreLFP(cycleMatrix(:,cycleNum));
			%figure; plot(startToEndCycle,'r*')
			%figure; plot(rawStartToEndCyclePadded((resamplePadLength+1):(end-resamplePadLength)),'r*')
			%figure; plot(resampledCycle(2:(end-1)),'b*')	
			%figure; plot(cycleMatrix(:,cycleNum),'b*')	
			%fds
		    else
			cycleMatrix(:,cycleNum)=decLFP(cycleStartIdx:cycleEndIdx)/cycleMax;
		    end
			totalNumCyclesSaved=totalNumCyclesSaved+1;
			%totalNumCyclesSaved/numCycles	
		%sprintf('%.0f%',totalNumCyclesSaved/numCycles*100)
			actualDurations(totalNumCyclesSaved)=cycleDuration;
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
%save(sprintf('cycleMatrixBroadbandZscore-MG49-%s-%s-%0.1f-ms-Window.mat',freqBandName,zscoreStr,halfCycleWindow*2000),'cycleMatrix','-v7.3')
save(sprintf('cycleMatrixBroadbandZscore-MG49-%s-%s-DurationNormalized.mat',freqBandName,zscoreStr),'actualDurations','cycleMatrix','-v7.3')

toc

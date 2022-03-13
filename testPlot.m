filePath=saveFileName;
data=load(filePath);
close all
trackLength=120; %cm
positions=data.spaceTimePhaseInfo.spaceTimePhasePerCycle{:,1}; %normalized


times=data.spaceTimePhaseInfo.spaceTimePhasePerCycle{:,2}; %seconds?

speeds=data.spaceTimePhaseInfo.spaceTimePhasePerCycle{:,5};

timeSpentAtEachPosition=data.spikePerCycleInfo.timeSpentAtEachPosition;

plot(linspace(0,1,length(timeSpentAtEachPosition)),timeSpentAtEachPosition,'k-')


phases=data.spaceTimePhaseInfo.spaceTimePhasePerCycle{:,3};

cycleDists=data.spaceTimePhaseInfo.spaceTimePhasePerCycle{:,4};

cycleDurations=data.spikePerCycleInfo.durationPerCycleSec;


topSpeed=prctile(data.spikePerCycleInfo.speedPerCycle,95);
trackMaxTime=trackLength/topSpeed; %sec

figure
subplot(2,2,3)

linearPos=[];
linearTime=[];
phasesInRadians=[];

for pi=1:length(positions)
	currPhases=phases(pi,:);
	currPhases(isnan(currPhases))=[];
	
	currPositions=repelem(positions(pi),length(currPhases));
	currTimes=repelem(times(pi),length(currPhases));
        for si=1:length(currPhases)
		currPositions(si)=currPositions(si)+currPhases(si)/360*cycleDists(pi);
		currTimes(si)=currTimes(si)+currPhases(si)/360*cycleDurations(pi);
		linearPos=[linearPos currPositions(si)];
		linearTime=[linearTime currTimes(si)];
		phasesInRadians=[phasesInRadians currPhases(si)*pi/180];
	end
	
	plot(currPositions*trackLength,currPhases,'k.')
	hold on
end
xlabel('Distance along track (cm)')
ylabel('Phase in theta cycle (deg)')
xlim([0 trackLength])
ylim([0 360])


cdfSpaceBinSize=1; %cm
%cdfTimeBinSize=1; %sec
numCDFbins=trackLength/cdfSpaceBinSize;
cdfSpikesVsPos=NaN(numCDFbins,1)
cdfSpikesVsSpaceTime=NaN(numCDFbins,1)

pdfSpikesVsSpaceTime=NaN(numCDFbins,numCDFbins);

positionsCDF=linspace(1,0,numCDFbins);
timesCDF=linspace(1,0,numCDFbins);

for ci=1:numCDFbins
	currPos=positionsCDF(ci);
	currTime=timesCDF(ci)*trackMaxTime;

	allSpikePositionsBeyondCurrPos=linearPos(linearPos>currPos);
	allSpikePositionsBeyondCurrPosAndTime=linearPos(linearPos>currPos & linearTime>currTime);
	cdfSpikesVsPos(ci)=length(allSpikePositionsBeyondCurrPos)/length(phasesInRadians);
	%cdfSpikesVsSpaceTime(ci)=length(allSpikePositionsBeyondCurrPosAndTime)/length(phasesInRadians);

end

cdfSpikesVsSpaceTime=spaceTimePhaseInfo.spaceTimeCDF(:);
cdfSpikesVsSpaceTime=flipud(cdfSpikesVsSpaceTime);

plot(positionsCDF*trackLength,cdfSpikesVsPos*360,'r-','LineWidth',3)
hold on
%plot(positionsCDF*trackLength,cdfSpikesVsSpaceTime*360,'b-','LineWidth',3)


subplot(2,2,1)
%histogram(positions*trackLength,30)
goodCycles=data.spaceTimePhaseInfo.goodCyclesForCell;

%if(di==1)
%goodPositions=(1-data.spikePerCycleInfo.normPosPerCycle(goodCycles))*trackLength;
%else
goodPositions=(data.spikePerCycleInfo.normPosPerCycle(goodCycles))*trackLength;
%end
goodSpikeRates=data.spikePerCycleInfo.numSpikesPerCycle(goodCycles);

spaceBinSize=3;

numBins=(trackLength/spaceBinSize);
centerPos=NaN(numBins,1);
currSpikeRate=zeros(numBins,1);
for bi=1:numBins
	currStartPos=(bi-1)*spaceBinSize;
	centerPos(bi)=currStartPos+spaceBinSize/2;
	currEndPos=(bi)*spaceBinSize;
	currBinIdxes=goodPositions>=currStartPos & goodPositions<currEndPos;
	if(~isempty(currBinIdxes))
		currSpikeRate(bi)=nanmean(goodSpikeRates(currBinIdxes));
	end
end

%plot(goodPositions, goodSpikeRates,'-k')
%plot(centerPos,currSpikeRate,'-k')
if(di==1)
	plot(flipud(spaceTimePhaseInfo.placeField(:)),'k-')
	fieldStart=trackLength-spaceTimePhaseInfo.startFieldCm;
	fieldEnd=trackLength-spaceTimePhaseInfo.endFieldCm;
else
	plot(spaceTimePhaseInfo.placeField(:),'k-')
	fieldStart=spaceTimePhaseInfo.startFieldCm;
	fieldEnd=spaceTimePhaseInfo.endFieldCm;
end
xlim([0 trackLength])
xlabel('Distance along track (cm)')
ylabel('Firing rate (Hz)')
hold on
plot([fieldStart fieldStart],ylim,'b')
plot([fieldEnd fieldEnd],ylim,'r')

subplot(2,2,2)
histogram(data.spaceTimePhaseInfo.spaceTimePhasePerCycle{:,4}*trackLength)
xlabel('Distance travelled during theta cycle (cm)')
ylabel('Number of theta cycles')
xlim([0 25])

size(data.spaceTimePhaseInfo.spaceTimePhasePerCycle{:,4})

size(data.spikePerCycleInfo.speedPerCycle(goodCycles))
subplot(2,2,4)
plot(positions*trackLength,speeds,'.k')

%histogram(data.spikePerCycleInfo.speedPerCycle(goodCycles))
%xlabel('Speed (cm/sec)')
%%%%ylabel('Number of theta cycles')
xlabel('Distance along track (cm)')
ylabel('Speed in field (cm/sec)')
xlim([0 trackLength])
maxFig
setFigFontTo(18)
uberTitle(removeUnderscores(fileRootName))

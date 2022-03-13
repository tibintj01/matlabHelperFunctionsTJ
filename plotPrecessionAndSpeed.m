
data=currData;
spaceTimePhaseInfo=currData.spaceTimePhaseInfo;
close all
trackLength=120; %cm
positions=data.spaceTimePhaseInfo.spaceTimePhasePerCycle{:,1}; %normalized
times=data.spaceTimePhaseInfo.spaceTimePhasePerCycle{:,2}; %seconds?
speeds=data.spaceTimePhaseInfo.spaceTimePhasePerCycle{:,5};

timeSpentAtEachPosition=data.spikePerCycleInfo.timeSpentAtEachPosition;


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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %adjust phase based on manual offset (to correct for wrong theta
    %reference
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	currPhases= mod(currPhases-(manualPrecessionStartPhaseDeg-360),360);
	plot(currPositions*trackLength,currPhases,'k.')
	hold on
end
xlabel('Distance along track (cm)')
ylabel('Phase in theta cycle (deg)')
xlim([currStartCm currEndCm])
ylim([0 360])
axis square

title('Spike Phase vs Distance')

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

%plot(positionsCDF*trackLength,cdfSpikesVsPos*360,'r-','LineWidth',3)
hold on
%plot(positionsCDF*trackLength,cdfSpikesVsSpaceTime*360,'b-','LineWidth',3)


subplot(2,2,2)
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
plotFiringRateField



subplot(2,2,4)
plot(positions*trackLength,speeds,'ok','MarkerSize',3)
hold on
plot([currStartCm currStartCm],[0 40],'b','LineWidth',1)
plot([currEndCm currEndCm],[0 40],'r','LineWidth',1)
title('Speed vs Distance')

%histogram(data.spikePerCycleInfo.speedPerCycle(goodCycles))
%xlabel('Speed (cm/sec)')
%%%%ylabel('Number of theta cycles')
xlabel('Distance along track (cm)')
ylabel('Speed in field (cm/sec)')
xlim([0 trackLength])
ylim([0 40])
maxFig
setFigFontTo(18)
uberTitle(removeUnderscores(saveRootName))

subplot(2,2,1)
plotFiringRateField
xlim([currStartCm currEndCm])
axis square
setFigFontTo(18)


linearPos=[];
linearTime=[];
phasesInRadians=[];
useManualFieldBoundaries=0;
positions=data.spaceTimePhaseInfo.spaceTimePhasePerCycle{:,1}; %normalized

phases=data.spaceTimePhaseInfo.spaceTimePhasePerCycle{:,3};
times=data.spaceTimePhaseInfo.spaceTimePhasePerCycle{:,2}; %seconds?
trackLength=120; %cm
cycleDists=data.spaceTimePhaseInfo.spaceTimePhasePerCycle{:,4};
cycleDurations=data.spikePerCycleInfo.durationPerCycleSec;
topSpeed=prctile(data.spikePerCycleInfo.speedPerCycle,95);
trackMaxTime=trackLength/topSpeed; %sec

manualPrecessionStartPhaseDeg=data.manualPrecessionStartPhaseDeg;
fieldStart=data.manualFieldStartCm;
fieldEnd=data.manualFieldEndCm;

allSpeedChanges=[];
allPhases=[];
allDists=[];
allTimes=[];
for posIdx=1:length(positions)
	currPhases=phases(posIdx,:);
	currPhases(isnan(currPhases))=[];
	
	currPositions=repelem(positions(posIdx),length(currPhases));
	currTimes=repelem(times(posIdx),length(currPhases));
        for si=1:length(currPhases)
		currPositions(si)=currPositions(si)+currPhases(si)/360*cycleDists(posIdx);
		currTimes(si)=currTimes(si)+currPhases(si)/360*cycleDurations(posIdx);
		linearPos=[linearPos currPositions(si)];
		linearTime=[linearTime currTimes(si)];
		phasesInRadians=[phasesInRadians currPhases(si)*posIdx/180];
     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %adjust phase based on manual offset (to correct for wrong theta
    %reference
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	currPhases= mod(currPhases-(manualPrecessionStartPhaseDeg-360),360);
    currDists=currPositions*trackLength;
    
    if(useManualFieldBoundaries)
        currInFieldIdxes=currDists>=fieldStart & currDists <=fieldEnd;
        currDists=currDists(currInFieldIdxes);
        currPhases=currPhases(currInFieldIdxes);
    end
  
   plot(currDists,currPhases,'k.','MarkerSize',10)

	hold on
    allPhases=[allPhases;currPhases(:)];
    allDists=[allDists; currDists(:)];
    allTimes=[allTimes; currTimes(:)];
end

    [circLinRhoDist,circLinPDist,circLinSlopeDist,circLinOffsetDist]=kempter_lincirc(allDists,ang2rad(allPhases));
xlabel('Distance along track (cm)')
ylabel('Phase in theta cycle (deg)')
    ylim([0 360])
    circLinSlopeDist=circLinSlopeDist*360;
    

    title({'Phase vs distance',sprintf('circlinear: rho = %.4f, slope= %.2f, p=%.2f',circLinRhoDist,circLinSlopeDist,circLinPDist)})

%title('SposIdxke Phase vs Distance')

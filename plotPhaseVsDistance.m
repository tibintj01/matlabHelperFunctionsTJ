linearPos=[];
linearTime=[];
phasesInRadians=[];
positions=data.spaceTimePhaseInfo.spaceTimePhasePerCycle{:,1}; %normalized

phases=data.spaceTimePhaseInfo.spaceTimePhasePerCycle{:,3};
times=data.spaceTimePhaseInfo.spaceTimePhasePerCycle{:,2}; %seconds?
trackLength=120; %cm
cycleDists=data.spaceTimePhaseInfo.spaceTimePhasePerCycle{:,4};
cycleDurations=data.spikePerCycleInfo.durationPerCycleSec;
topSpeed=prctile(data.spikePerCycleInfo.speedPerCycle,95);
trackMaxTime=trackLength/topSpeed; %sec

manualPrecessionStartPhaseDeg=data.manualPrecessionStartPhaseDeg;

allSpeedChanges=[];
allPhases=[];
allDists=[];
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
    currInFieldIdxes=currDists>=fieldStart & currDists <=fieldEnd;
    currDists=currDists(currInFieldIdxes);
    currPhases=currPhases(currInFieldIdxes);
    
    if(convertDistanceToSpeedChange)
        %currSpeedChanges=-interp1(speedChangeDists,speedChanges,currDists,'linear','extrap');
         currSpeedChanges=-interp1(speedChangeDists,speedChanges,currDists,'linear');
        
        
        %plot(currSpeedChanges,currPhases,'b.','MarkerSize',10)
       
       
        allSpeedChanges=[allSpeedChanges; currSpeedChanges(:)];
        
        
        
    else
        %plot(currDists,currPhases,'k.','MarkerSize',10)
    end
	hold on
    allPhases=[allPhases;currPhases(:)];
    allDists=[allDists; currDists(:)];
end




if(convertDistanceToSpeedChange)
        allReasonableSpeedIdxes=(abs(allSpeedChanges(:))<maxReasonableSpeedChange);
    allReasonableSpeedChanges=allSpeedChanges(allReasonableSpeedIdxes);
    allReasonablePhases=allPhases(allReasonableSpeedIdxes);
    [circLinRhoSpeed,circLinPSpeed,circLinSlopeSpeed,circLinOffsetSpeed]=kempter_lincirc(allReasonableSpeedChanges,ang2rad(allReasonablePhases));
    if(showPlots)
     xlabel('-Speed change along track (cm/s per cm)')
     plot(allReasonableSpeedChanges,allReasonablePhases,'b.','MarkerSize',10)
    end
else
    allReasonableSpeedDistanceIdxes=allDists>maxReasonableSpeedCorrespondingDist & allDists<trackLength-maxReasonableSpeedCorrespondingDist;
    allReasonableDists=allDists(allReasonableSpeedDistanceIdxes);
    allReasonablePhases=allPhases(allReasonableSpeedDistanceIdxes);
    [circLinRhoDist,circLinPDist,circLinSlopeDist,circLinOffsetDist]=kempter_lincirc(allReasonableDists,ang2rad(allReasonablePhases));
    if(showPlots)
        plot(allReasonableDists,allReasonablePhases,'k.','MarkerSize',10)
        xlabel('Distance along track (cm)')
    end
end
try
    circLinSlopeSpeed=circLinSlopeSpeed*360;
end
try
    circLinSlopeDist=circLinSlopeDist*360;
end
if(showPlots)
    ylabel('Phase in theta cycle (deg)')
    ylim([0 360])
    axis square

end

%title('SposIdxke Phase vs Distance')

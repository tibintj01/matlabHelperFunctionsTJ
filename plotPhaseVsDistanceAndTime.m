linearPos=[];
linearTime=[];
phasesInRadians=[];
showPlots=0;
useManualFieldBoundaries=0;
useManualFieldBoundaries=1;
positions=data.spaceTimePhaseInfo.spaceTimePhasePerCycle{:,1}; %normalized
speeds=data.spaceTimePhaseInfo.spaceTimePhasePerCycle{:,5}; 


phases=data.spaceTimePhaseInfo.spaceTimePhasePerCycle{:,3};
times=data.spaceTimePhaseInfo.spaceTimePhasePerCycle{:,2}; %seconds?
trackLength=120; %cm



cycleDists=data.spaceTimePhaseInfo.spaceTimePhasePerCycle{:,4};
cycleDurations=data.spikePerCycleInfo.durationPerCycleSec;
topSpeed=prctile(data.spikePerCycleInfo.speedPerCycle,95);
trackMaxTime=trackLength/topSpeed; %sec

try
manualPrecessionStartPhaseDeg=data.manualPrecessionStartPhaseDeg;
fieldStart=data.manualFieldStartCm;
fieldEnd=data.manualFieldEndCm;
findBestShift=0;
catch
    findBestShift=1;
    fieldStart=0;
   fieldEnd=120;
end

allPhases=[];
allDists=[];
allTimes=[];
allSpeeds=[];
for posIdx=1:length(positions)
	currPhases=phases(posIdx,:);
	currPhases(isnan(currPhases))=[];
	
	currPositions=repelem(positions(posIdx),length(currPhases));
	currTimes=repelem(times(posIdx),length(currPhases));
    currSpeeds=repelem(speeds(posIdx),length(currPhases));
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
        currTimes=currTimes(currInFieldIdxes);
        currSpeeds=currSpeeds(currInFieldIdxes);
    end
    
  % plot(currDists,currPhases,'k.','MarkerSize',10)

	%hold on
    allPhases=[allPhases;currPhases(:)];
    allDists=[allDists; currDists(:)];
    allTimes=[allTimes; currTimes(:)];
    allSpeeds=[allSpeeds; currSpeeds(:)];
end

if(findBestShift)
 [bestPhaseShift]=getBestLinearCorrShift(allPhases,allDists);
    phases=mod(phases+bestPhaseShift,360);
    
end
   %for combining fields
   minDist=min(allDists);
   maxDist=max(allDists);

   %zscoreTimes=zscoreLFP(allTimes);
   %outlierTimeIdxes=abs(zscoreTimes)>3;
   %allTimes(outlierTimeIdxes)=NaN
   minSpeed=prctile(allSpeeds,0);
   maxTime=(maxDist)/minSpeed;
   
   
    %allTimes(allTimes>maxTime)=NaN;
   minTime=min(allTimes);
   maxTime=max(allTimes);

    
     allTimes=(allDists-fieldStart)./allSpeeds;
   allDists=scaledata(allDists,0,1);
   %allTimes=scaledata(allTimes,0,1);
    %allTimes=zscoreLFP(allTimes);
     %allTimes=allTimes-min(allTimes);
    
     
     if(showPlots)

    [circLinRhoDist,circLinPDist,circLinSlopeDist,circLinOffsetDist]=kempter_lincirc(allDists,ang2rad(allPhases));
    [circLinRhoTime,circLinPTime,circLinSlopeTime,circLinOffsetTime]=kempter_lincirc(allTimes,ang2rad(allPhases));
    [circLinRhoDT,circLinPDT,circLinSlopeDT,circLinOffsetDT]=kempter_lincirc(allTimes.*allDists,ang2rad(allPhases));
    subplot(3,1,1)
 plot(allDists,allPhases,'k.','MarkerSize',10)
    %xlabel('Distance along track (cm)')
    xlabel('Fraction distance in field')
ylabel('Phase in theta cycle (deg)')
    ylim([0 360])
    circLinSlopeDist=circLinSlopeDist*360;
    title({sprintf('Phase vs distance (field: %.1f to %.1f cm)',minDist,maxDist),sprintf('circlinear: rho = %.4f, slope= %.2f, p=%.2f',circLinRhoDist,circLinSlopeDist,circLinPDist)})

    
    subplot(3,1,2)
     plot(allTimes,allPhases,'k.','MarkerSize',10)
     xlim([0 5])
    %xlabel('Time elapsed in lap (sec)')
    xlabel('Zscore time elapsed in field (minus minimum)')
ylabel('Phase in theta cycle (deg)')
    ylim([0 360])
    circLinSlopeTime=circLinSlopeTime*360;
    title({sprintf('Phase vs time (field: %.2f to %.2f sec)',minTime,maxTime),sprintf('circlinear: rho = %.4f, slope= %.2f, p=%.2f',circLinRhoTime,circLinSlopeTime,circLinPTime)})
%title('SposIdxke Phase vs Timeance')

subplot(3,1,3)
     plot(allTimes.*allDists,allPhases,'k.','MarkerSize',10)
     xlim([0 5])
    %xlabel('Time elapsed in lap (sec)')
    xlabel('Fraction distance * zscore time')
ylabel('Phase in theta cycle (deg)')
    ylim([0 360])
    title({sprintf('Phase vs Dist*Time'),sprintf('circlinear: rho = %.4f, slope= %.2f, p=%.2f',circLinRhoDT,circLinSlopeDT,circLinPDT)})
    
        disp('here')
     end

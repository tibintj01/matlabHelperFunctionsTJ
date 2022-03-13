function [spaceTimePhaseInfo] = getSpaceTimePhaseInfo(spikePerCycleInfo,dirFlag)
    if(dirFlag<1)
        dirFlag=-1;
        dirID=1;
    else
        dirFlag=1;
        dirID=2;
    end
    
     normPosFieldExpectedStart=spikePerCycleInfo.peakFieldStartPosNorm(dirID);
     normPosFieldExpectedEnd=spikePerCycleInfo.peakFieldEndPosNorm(dirID);
     if(isnan(normPosFieldExpectedEnd) && isnan(normPosFieldExpectedStart))
        normPosFieldExpectedStart=spikePerCycleInfo.peakFieldStartPosNorm(2);
        normPosFieldExpectedEnd=spikePerCycleInfo.peakFieldEndPosNorm(2);
     end
    cellIDstr=spikePerCycleInfo.cellIDstr;

    getSpaceTimeProps=0;
    
    %nanOutsideField=1;
    nanOutsideField=0;
    combineFields=1;
    %savePlots=1;
    savePlots=0;
    shift180=1;
    shift180=0;

    adjustedTimes=spikePerCycleInfo.elapsedTimeInLapPerCycle;
    %maxTime=prctile(adjustedTimes,85);
    %maxTime=20;
    maxTime=10;
    adjustedTimes(adjustedTimes<0)=0;
    %adjustedTimes(adjustedTimes>maxTime)=maxTime;

    %dir Flag is 1 or -1 
    fwdLapIdxes=spikePerCycleInfo.lapDirectionPerCycle==dirFlag;

    %minSpeedForKeepingSpike=10;
    minSpeedForKeepingSpike=5;
    %minSpeedForKeepingSpike=1;
    %minSpeedForKeepingSpike=0;
    %minThetaZampForKeepingSpike=0;
    minThetaZampForKeepingSpike=2.5;
    %minThetaZampForKeepingSpike=2;
    %minThetaZampForKeepingSpike=1.5;
    %minThetaZampForKeepingSpike=-Inf;

    movingIdxes=spikePerCycleInfo.speedPerCycle>minSpeedForKeepingSpike;
    highThetaIdxes=spikePerCycleInfo.thetaAmpZPerCycle >minThetaZampForKeepingSpike;   


    goodCyclesForCell=fwdLapIdxes(:) & movingIdxes(:) & highThetaIdxes(:);     
    %inFieldCycleIdxes=spikePerCycleInfo.isInFieldPerCycle(~fwdLapIdxes,dirID);
    %inFieldCycleIdxes=spikePerCycleInfo.isInFieldPerCycle(fwdLapIdxes,dirID);
    %inFieldCycleIdxes=spikePerCycleInfo.isInFieldPerCycle(goodCyclesForCell,dirID);
    inFieldCycleIdxes=spikePerCycleInfo.isInFieldPerCycle(:,dirID);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %get firing rate vs position curve
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    goodCyclePositions=spikePerCycleInfo.normPosPerCycle(goodCyclesForCell);
    goodCycleFiringRates=spikePerCycleInfo.firingRatePerCycle(goodCyclesForCell);


    
    if(nanOutsideField)
    	goodCyclesForCell=goodCyclesForCell & inFieldCycleIdxes(:);
    end

%figure; plot(spikePerCycleInfo.normPosPerCycle(goodCyclesForCell),spikePerCycleInfo.firingRatePerCycle(goodCyclesForCell),'.'); ylabel('firing rate(Hz)'); xlabel('position')


    firingRateVsPosWindow=10; %cm
    windowShiftSize=1; %cm
    trackLength=120; %cm
    numPointsPlaceField=trackLength/windowShiftSize;
    placeField=NaN(numPointsPlaceField,1);
    
    startFieldCm=spikePerCycleInfo.peakFieldStartPosNorm(dirID) * trackLength;
    endFieldCm=spikePerCycleInfo.peakFieldEndPosNorm(dirID) * trackLength;
    
    for wi=1:windowShiftSize:numPointsPlaceField
	currWindPosStart=(wi-1);
	currWindPosEnd=currWindPosStart+firingRateVsPosWindow;

        cyclesInCurrWindow=goodCyclePositions>currWindPosStart/trackLength & goodCyclePositions<currWindPosEnd/trackLength;

	firingRatesInCurrWindow=goodCycleFiringRates(cyclesInCurrWindow);

	placeField(wi)=nanmean(firingRatesInCurrWindow);
    end

   %figure; plot(1:numPointsPlaceField,placeField,'k-')
	%xlabel('Position (cm)')
	%ylabel('Firing rate (Hz)')
    
    spaceTimeCDF=NaN(size(placeField));
    for pi=length(placeField):-1:1
        ratesAboveCurrentPos=placeField(pi:end);
        spaceTimeCDF(pi)=nansum(ratesAboveCurrentPos);
    end
    
    spaceTimeCDF=spaceTimeCDF/max(spaceTimeCDF);
    %figure; plot(1:numPointsPlaceField,spaceTimeCDF)

    
    cycleDistances=spikePerCycleInfo.distTravelledPerCycleNormTrack(goodCyclesForCell);
    
    speedsInField=spikePerCycleInfo.speedPerCycle(goodCyclesForCell);
    %spaceTimePhasePerCycle=[spikePerCycleInfo.posPerCycle(~fwdLapIdxes),...
    %spaceTimePhasePerCycle=[spikePerCycleInfo.normPosPerCycle(~fwdLapIdxes),...
    %    adjustedTimes(~fwdLapIdxes), ...
    %    spikePerCycleInfo.firstSpikePhasePerCycle(~fwdLapIdxes)];
    %spaceTimePhasePerCycle=[spikePerCycleInfo.normPosPerCycle(fwdLapIdxes),...
    %    adjustedTimes(fwdLapIdxes), ...
    %    spikePerCycleInfo.firstSpikePhasePerCycle(fwdLapIdxes)];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %flip space for R to L direction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    positions=spikePerCycleInfo.normPosPerCycle(goodCyclesForCell);
    if(dirFlag==-1)
        positions=max(positions)-positions;
    end

    elapsedTimes=adjustedTimes(goodCyclesForCell);   

    spikePhasesInThisCycle=[spikePerCycleInfo.allSpikesPerCycle(goodCyclesForCell,:)];
 
    spaceTimePhasePerCycle=table(positions,...
        elapsedTimes, ...
        spikePhasesInThisCycle,...
		cycleDistances,speedsInField);
        %spikePerCycleInfo.firstSpikePhasePerCycle(goodCyclesForCell)];


     avgSpeed=spikePerCycleInfo.meanSpeed;
     timeFieldExpectedStart=normPosFieldExpectedStart*avgSpeed;
     timeFieldExpectedEnd=normPosFieldExpectedEnd*avgSpeed;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %DISCRETIZE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    discretizeSpaceAndTime=0;
    if(discretizeSpaceAndTime)
	    numBinsSpace=50;
	    numBinsTime=15;
	     %numBinsTime=50;
	     
	    spaceTimePhasePerCycle(:,1)=round(scaledata(spaceTimePhasePerCycle(:,1),0,numBinsSpace));
	    spaceTimePhasePerCycle(:,2)=round(scaledata(spaceTimePhasePerCycle(:,2),0,numBinsTime));
    else
	    %spaceTimePhasePerCycle(:,1)=scaledata(spaceTimePhasePerCycle(:,1),0,1);
	    %spaceTimePhasePerCycle(:,2)=scaledata(spaceTimePhasePerCycle(:,2),0,1);
	    spaceTimePhasePerCycle(:,1)=table(scaledata(spaceTimePhasePerCycle{:,1},0,1));
	    %spaceTimePhasePerCycle(:,2)=table(scaledata(spaceTimePhasePerCycle{:,2},0,1)); %
    end
    
    if(nanOutsideField)
        %spaceTimePhasePerCycle(~inFieldCycleIdxes,:)=[];
    end
    disp('testing')
    %spaceTimePhasePerCycle(any(isnan(spaceTimePhasePerCycle),2),:)=[];
    spaceTimePhasePerCycle(any(isnan(spaceTimePhasePerCycle{:,1:2}),2),:)=[];

    badRows=[];
    for ri=1:size(spaceTimePhasePerCycle,1)
        currSpikeData=spaceTimePhasePerCycle{ri,3};
        if(isnan(max(currSpikeData)))
            badRows=[badRows; ri];
        end
    end
    spaceTimePhasePerCycle(badRows,:)=[];
    
    if(getSpaceTimeProps)
	    %zs=threeDPtsToSurface(round(spaceTimePhasePerCycle'));
	    zs=threeDPtsToSurface((spaceTimePhasePerCycle'));
	    avgPhasePerSpaceTimeBin=NaN(size(zs));
	    for row=1:size(zs,1)
		for col=1:size(zs,2)
		    %if(~isempty(zs{row,col}))
		    if(length(zs{row,col})>0)
		     avgPhasePerSpaceTimeBin(row,col)=circMeanDeg(zs{row,col});
		    end
		end
	    end

	    %spaceBins=unique(round(spaceTimePhasePerCycle(:,1)))/numBinsSpace;
	    spaceBins=unique(round(spaceTimePhasePerCycle(:,1)));
	    spaceBins=scaledata(spaceBins,normPosFieldExpectedStart,normPosFieldExpectedEnd);
	    timeBins=unique(round(spaceTimePhasePerCycle(:,2)));
	   timeBins=scaledata(timeBins,timeFieldExpectedStart,timeFieldExpectedEnd);
	   if(shift180)
		    avgPhasePerSpaceTimeBin=avgPhasePerSpaceTimeBin-180;
	   end
		
	   avgPhasePerSpaceBin=circMeanDeg(avgPhasePerSpaceTimeBin,2);
	   avgPhasePerTimeBin=circMeanDeg(avgPhasePerSpaceTimeBin,1);
		if(savePlots) 
		  
			setTightSubplots_SpaceTime
			fH=figure
			subplot(2,2,2)
			omarPcolor(spaceBins,timeBins,avgPhasePerSpaceTimeBin',fH)
		%cmp=hsv
		%hsvInvert=[cmp(ceil(end/2):end,:); cmp(1:floor(end/2),:)];
			colormap(hsv)
		%colormap(hsvInvert)
			cb=colorbar('north')
			cb.Position=[0.05 0.25 0.4 0.05 ];
			xlabel('Pos (frac of track length)')
			ylabel('Time since leaving wall (sec)')
			ylabel(cb,'First spike theta phase (deg)')
		xlim([0 1])
		ylim([-0.5 maxTime])
		%caxis([0 360])
		%cb

			
			subplot(2,2,4)
			plot(spaceBins,avgPhasePerSpaceBin,'ko')
			hold on
			plot(spaceBins,avgPhasePerSpaceBin+360,'ko')
			%xlim([-Inf Inf])
		xlim([0 1])
			ylim([0 720])
			xlabel('Pos (frac of track length)')
			ylabel('First spike theta phase (deg)')

			
		
	      
			subplot(2,2,1)
			plot(avgPhasePerTimeBin,timeBins+0.5,'ko','MarkerSize',5) %binCenters
			hold on
			plot(avgPhasePerTimeBin+360,timeBins+0.5,'ko','MarkerSize',5)
			xlim([0 720])
			%ylim([-0.5 Inf])
		ylim([-0.5 maxTime])

			xlabel('First spike theta phase (deg)')
			ylabel('Time since leaving wall (sec)')
			setFigFontTo(16)
			uberTitle(sprintf('%s Direction %d, Spacetime Receptive Field',removeUnderscores(cellIDstr), dirFlag))
			maxFig
			saveas(gcf,sprintf('SpaceTimeField_%s_dir%d.tif',cellIDstr,dirFlag))
		end
	
		spaceTimePhaseInfo.avgPhasePerTimeBin=avgPhasePerTimeBin;
		spaceTimePhaseInfo.avgPhasePerSpaceTimeBin=avgPhasePerSpaceTimeBin;	
		spaceTimePhaseInfo.avgPhasePerSpaceTimeBinZ=zscoreLFP(avgPhasePerSpaceTimeBin);	
		spaceTimePhaseInfo.shift180=shift180;
    		spaceTimePhaseInfo.avgPhasePerSpaceBin=avgPhasePerSpaceBin;
		spaceTimePhaseInfo.spaceBins=spaceBins;
		spaceTimePhaseInfo.timeBins=timeBins;

	end	
	%save results into struct for further processing
	spaceTimePhaseInfo.spaceTimePhasePerCycle=spaceTimePhasePerCycle;

    spaceTimePhaseInfo.startFieldCm=startFieldCm;
    spaceTimePhaseInfo.endFieldCm=endFieldCm;
    
    spaceTimePhaseInfo.spaceTimeCDF=spaceTimeCDF;
    spaceTimePhaseInfo.placeField=placeField;
    spaceTimePhaseInfo.firingRateVsPosWindow=firingRateVsPosWindow;
    spaceTimePhaseInfo.windowShiftSize=windowShiftSize;
    
	spaceTimePhaseInfo.goodCyclesForCell=goodCyclesForCell;	
	spaceTimePhaseInfo.cellIDstr=cellIDstr;
    	spaceTimePhaseInfo.dirFlag=dirFlag;
    

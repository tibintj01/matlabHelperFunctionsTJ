    data=load(saveFileName);

    
	corrPThresh=0.01;

    minNumSpikes=100;
    skip=0;

    minRate=10;
    if(max(data.spikePerCycleInfo.peakFieldMaxRate)<minRate)
        skip=1;
        return
    end
    
	positions=data.spaceTimePhaseInfo.spaceTimePhasePerCycle{:,1};
	phases=data.spaceTimePhaseInfo.spaceTimePhasePerCycle{:,3};
	cycleDists=data.spaceTimePhaseInfo.spaceTimePhasePerCycle{:,4};
    
    if(isempty(positions) ||isempty(phases) ||isempty(cycleDists))
        skip=1;
        return
    end

	linearPos=[];
	phasesInRadians=[];

	for pi=1:length(positions)
		currPhases=phases(pi,:);
		currPhases(isnan(currPhases))=[];

		currPositions=repelem(positions(pi),length(currPhases));
		for si=1:length(currPhases)
			currPositions(si)=currPositions(si)+currPhases(si)/360*cycleDists(pi);
			linearPos=[linearPos currPositions(si)];
			phasesInRadians=[phasesInRadians currPhases(si)*pi/180];
		end
	end


    %if(isempty(linearPos) || isempty(phasesInRadians))
    if(length(linearPos)< minNumSpikes || length(phasesInRadians) < minNumSpikes)
        %skipCount=skipCount+1;
        skip=1;
    end

    [ rhoK,pK,sK,bK ] = kempter_lincirc( linearPos,phasesInRadians);

    numDensityBins=30;
    showInfoPlots=0;
    [IposPhase,excessJointSpikeDensity,Hpos,Hphase,rhoSpear,pSpear,dispXlims,dispYlims]=getDiscreteMutualInfo(linearPos,phasesInRadians,numDensityBins,numDensityBins,'position','phase (radians)',showInfoPlots);

    %if(pK>corrPThresh || isnan(pK))
    if(pSpear>corrPThresh || isnan(pK) || rhoSpear>0 || isnan(rhoSpear))
        %skipCount=skipCount+1;
        skip=1;
    end

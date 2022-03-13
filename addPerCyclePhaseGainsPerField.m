close all; clear all; clc
data=load('allUnitInFieldGoodSpeedSpikePhasesAndTimes.mat')


numFields=length(data.spikePhasesPerField);

phaseGainsPerCyclePerField=cell(numFields,1);
timeGainsPerCyclePerField=cell(numFields,1);
timeInFieldPerCyclePerField=cell(numFields,1);
startPhasePerCyclePerField=cell(numFields,1);

phaseKernelDistrBinCenters=1:1:360;
%gaussKernelWidth=45; %deg
gaussKernelWidth=15; %deg
%gaussKernelWidth=30; %deg

useCircMean=1;
fH=figure

for fi=1:numFields
    fi
    
    currFieldPhases=data.spikePhasesPerField{fi};
    currFieldSpikeTimesInExp=data.spikeTimesInExpPerField{fi};
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %ALREADY USES AUTOMATED TIME BOUND DEFINITION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    currFieldSpikeTimeFracInField=data.autoRenormedSpikeTimeFracInFieldPerField{fi};
    
    
    
    if(isempty(currFieldPhases) || isempty(currFieldSpikeTimesInExp))
        continue
    end
    
    currUnitInfo=load(data.unitInfoPathPerField{fi});
    currUnitCycleInfo=load(currUnitInfo.refChCycleInfoPath);
    
    ampZPerCycle=zscoreLFP(currUnitCycleInfo.ampPerCycle);
    
    ampZThresh=0.5;
    
    highAmpCycles=ampZPerCycle>=ampZThresh;
    
    thetaCycleMinTimes=currUnitCycleInfo.cycleMinTimes;
    thetaCycleMaxTimes=currUnitCycleInfo.cycleMaxTimes;
    
    numCycles=length(thetaCycleMinTimes);
    
    currFieldPhaseGainsPerCycle=NaN(numCycles,1);
    currFieldSummaryPhasePerCycle=NaN(numCycles,1);
    currFieldTimeInFieldFracPerCycle=NaN(numCycles,1);
    
    tic
    for ci=1:numCycles
        
        if(~highAmpCycles(ci))
            continue
        end
        currCycleMinTime=thetaCycleMinTimes(ci);
        
        if(currCycleMinTime<min(currFieldSpikeTimesInExp)-1)
            continue
        end
        
        [~,currCycleStartIdx]=min(abs(thetaCycleMaxTimes-currCycleMinTime));
        
        if(thetaCycleMaxTimes(currCycleStartIdx)>currCycleMinTime)
            currCycleStartIdx=currCycleStartIdx-1;
        end
        if(currCycleStartIdx==0 || currCycleStartIdx==length(thetaCycleMaxTimes))
            continue
        end
        
        currCycleStartTime=thetaCycleMaxTimes(currCycleStartIdx);
        currCycleEndTime=thetaCycleMaxTimes(currCycleStartIdx+1);
        
        currCycleSpikeTimeIdxes=currFieldSpikeTimesInExp>currCycleStartTime & currFieldSpikeTimesInExp<=currCycleEndTime;
        
        %currCycleSpikeTimeFracInCycle=(currFieldSpikeTimesInExp(currCycleSpikeTimeIdxes)-currCycleStartTime)/(currCycleEndTime-currCycleStartTime);
        
        currCycleSpikePhases=currFieldPhases(currCycleSpikeTimeIdxes);
        
        currCycleSpikeTimeInFieldFracs=currFieldSpikeTimeFracInField(currCycleSpikeTimeIdxes);
        
        if(isempty(currCycleSpikePhases))
            continue
        end

        if(length(currCycleSpikePhases)==1)
            currCycleSummaryPhase=currCycleSpikePhases;
        else
            if(useCircMean)
                currCycleSummaryPhase=nanmean(currCycleSpikePhases);
            else
             vfObservations=currCycleSpikePhases(:)*2*pi/360;
             vfPDFSamples=phaseKernelDistrBinCenters(:)*2*pi/360;
             vfDomain=[0 2*pi];
             fSigma=gaussKernelWidth*2*pi/360;
             [vfEstimate] = circ_ksdensity(vfObservations, vfPDFSamples, vfDomain, fSigma);

             currPdist=vfEstimate/sum(vfEstimate(:));
             %plot(phaseKernelDistrBinCenters,currPdist,'k-','LineWidth',3)

            [maxVal,maxID]=max(currPdist);
            
            maxIDs=find(abs(currPdist-maxVal)<0.0001);
            
            
            currCycleSummaryPhase=nanmean(phaseKernelDistrBinCenters(maxIDs));
             %currCycleSummaryPhase=circMeanDeg(phaseKernelDistrBinCenters(maxIDs));


            %currCycleSummaryPhase=phaseKernelDistrBinCenters(maxID);
            end
        end
        
        currFieldSummaryPhasePerCycle(ci)=currCycleSummaryPhase;
        
        currFieldTimeInFieldFracPerCycle(ci)=nanmean(currCycleSpikeTimeInFieldFracs);
        
    end
    toc
    %for ci=1:numCycles-1
    %{
   for ci=2:numCycles
        %currFieldPhaseGainsPerCycle(ci)=mod(currFieldSummaryPhasePerCycle(ci)-currFieldSummaryPhasePerCycle(ci-1)+180,360)-180;
        %currFieldTimeInFieldFracPerCycle(ci)=nanmean([currFieldTimeInFieldFracPerCycle(ci+1), currFieldTimeInFieldFracPerCycle(ci)]);
        currFieldTimeInFieldFracPerCycle(ci)=nanmean([currFieldTimeInFieldFracPerCycle(ci), currFieldTimeInFieldFracPerCycle(ci-1)]);

        currFieldPhaseGainsPerCycle(ci)=angdiffDeg( [ currFieldSummaryPhasePerCycle(ci-1) currFieldSummaryPhasePerCycle(ci) ]); 
   end
    %}
    
   %currFieldPhaseGainsPerCycle=angdiffDeg(currFieldSummaryPhasePerCycle);
   currFieldPhaseGainsPerCycle=diff(currFieldSummaryPhasePerCycle);

   currFieldTimeGainsPerCycle=diff(currFieldTimeInFieldFracPerCycle);
   
   currFieldPhaseGainsPerCycle=[currFieldPhaseGainsPerCycle(:); NaN];
   currFieldTimeGainsPerCycle=[currFieldTimeGainsPerCycle(:); NaN];
   
    notNaNIdxes=~isnan(currFieldPhaseGainsPerCycle) & ~isnan(currFieldTimeInFieldFracPerCycle);
    
    
    currFieldSummaryPhasePerCycle(~notNaNIdxes)=[];
    currFieldPhaseGainsPerCycle(~notNaNIdxes)=[];
    currFieldTimeInFieldFracPerCycle(~notNaNIdxes)=[];
    currFieldTimeGainsPerCycle(~notNaNIdxes)=[];
    
    
    %{
    plot(currFieldTimeInFieldFracPerCycle,1./currFieldPhaseGainsPerCycle,'k.','MarkerSize',10)
    hold on
    xlim([0 1])
    ylim([-0.1 0])
    drawnow
    
    figure(fH); 
    plot(currFieldSummaryPhasePerCycle,currFieldPhaseGainsPerCycle,'ko')
    hold on
    %}
    %plot(currFieldSummaryPhasePerCycle) 
    
    
    %assert(length(currFieldSummaryPhasePerCycle)==length(currFieldPhaseGainsPerCycle))
    
    timeGainsPerCyclePerField{fi}=currFieldTimeGainsPerCycle;
    startPhasePerCyclePerField{fi}=currFieldSummaryPhasePerCycle;
    phaseGainsPerCyclePerField{fi}=currFieldPhaseGainsPerCycle;
    timeInFieldPerCyclePerField{fi}=currFieldTimeInFieldFracPerCycle;
end

save('allUnitInFieldGoodSpeedSpikePhasesAndTimes.mat','timeGainsPerCyclePerField','startPhasePerCyclePerField','phaseGainsPerCyclePerField','timeInFieldPerCyclePerField', '-append')

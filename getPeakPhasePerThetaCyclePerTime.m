function [localSpikePhasePerTime]=getPeakPhasePerThetaCyclePerTime(gaussTimeAxis,localGaussSpikeRatePerTime,cycleStartTimes,newTimeAxis)
%UNTITLED4 Summary of this function goes here

    peakPhasePerCycle=NaN(size(cycleStartTimes));
    for ci=1:(length(cycleStartTimes))
       	if(mod(ci,1000)==0)
		 ci/length(cycleStartTimes)
	end
        currCycleStartTime=cycleStartTimes(ci);
        if(ci<length(cycleStartTimes))
            currCycleEndTime=cycleStartTimes(ci+1);
        else
            currCycleEndTime=min(newTimeAxis(end),cycleStartTimes(ci)+0.3);
        end
        
        currCycleGaussIdxes=gaussTimeAxis<currCycleEndTime & gaussTimeAxis>=currCycleStartTime;
        
        currCycleRates=localGaussSpikeRatePerTime(currCycleGaussIdxes);
        if(isempty(currCycleRates))
            continue
        end 
        
        [pkRates,pkIdxesInCycle]=findpeaks(currCycleRates);
        [maxRate,maxRateIdxIdx]=max(pkRates);
        
      
        %pick the representative phase as the highest local peak
        %if there are no peaks pick the center of mass as phase
        %if there is no rate keep as NaN
        if(isempty(pkRates))
            if(maxRate==0)
              continue
            end
            
            rateCOMinCycle=getCOM(currCycleRates);
            comRatioInCycle=rateCOMinCycle/length(currCycleRates);
            
            peakPhasePerCycle(ci)=comRatioInCycle*360;
            continue
        end
        
   
        %if there is a tie pick the circ average amongst them as phase
        if(length(find(pkRates==maxRate))==1)
            peakRateIdxInCycle=pkIdxesInCycle(1);
            peakIdxRatioInCycle=peakRateIdxInCycle/length(currCycleRates);
            peakPhasePerCycle(ci)=peakIdxRatioInCycle*360;
            
        elseif(length(find(pkRates==maxRate))>1)
            peakIdxRatiosInCycle=pkIdxesInCycle/length(currCycleRates);
            
            peakPhases=peakIdxRatiosInCycle*360;
            
            peakPhasePerCycle(ci)=circMeanDeg(peakPhases);
        end
        
        peakPhasePerCycle(ci)=mod(peakPhasePerCycle(ci)+180,360); %to match previous convention
    end

    localSpikePhasePerTime=NaN(size(newTimeAxis));
    for ti=1:length(newTimeAxis)
        currTime=newTimeAxis(ti);
		[~,closestCycleStartID]=min(abs(cycleStartTimes-currTime));

		closestCycleStartTime=cycleStartTimes(closestCycleStartID);

		%closest theta cycle start is before or at current time
		if(closestCycleStartTime<=currTime)
            closestCycleStartID=closestCycleStartID;
		%closest theta cycle start is after current time
		else
			closestCycleStartID=max(1,closestCycleStartID-1);
        end
        
        localSpikePhasePerTime(ti)=peakPhasePerCycle(closestCycleStartID);
    end



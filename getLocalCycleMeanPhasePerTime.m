function [localAvgPhasePerTime,localAvgTimingPerTime] =getLocalCycleMeanPhasePerTime(allSpikeTimesAndPhases,timeBinCenters,cycleStartTimes)
%gets local running circular avg phase matching timebins in timeBinCenters
%based on cycle times - with threshold on consistency

        numTimeBins=length(timeBinCenters);

        approxTimeStep=median(diff(timeBinCenters));
        spikeTimes=allSpikeTimesAndPhases(:,1);
        spikePhases=allSpikeTimesAndPhases(:,2);

        localAvgPhasePerTime=NaN(size(timeBinCenters));
        localAvgTimingPerTime=NaN(size(timeBinCenters));
        
	for ti=1:numTimeBins
		currTime=timeBinCenters(ti);
		[~,closestCycleStartID]=min(abs(cycleStartTimes-currTime));

		closestCycleStartTime=cycleStartTimes(closestCycleStartID);

		%closest theta cycle start is before or at current time
		if(closestCycleStartTime<=currTime)
			localWindStartTime=closestCycleStartTime;
			if(closestCycleStartID<length(cycleStartTimes))
				localWindEndTime=cycleStartTimes(closestCycleStartID+1);
			else
				localWindEndTime=max(timeBinCenters)+approxTimeStep/2;
			end

		%closest theta cycle start is after current time
		else
			localWindEndTime=closestCycleStartTime;
                        if(closestCycleStartID>1)
                                localWindStartTime=cycleStartTimes(closestCycleStartID-1);
                        else
                                localWindStartTime=0;
                        end
		end

                
                localSpikePhases=spikePhases(spikeTimes>=localWindStartTime & spikeTimes <localWindEndTime);
                localSpikeTimes=spikeTimes(spikeTimes>=localWindStartTime & spikeTimes <localWindEndTime);
                
                

                localSpikePhases(isnan(localSpikePhases))=[];
                localSpikeTimes(isnan(localSpikeTimes))=[];
                
                
                
                if(~isempty(localSpikePhases))
                    %{
                    if(length(localSpikePhases)>3 && circ_r(ang2rad(localSpikePhases))<0.5)
                        figure;
                        polarhistogram(ang2rad(localSpikePhases))
                        disp(circMeanDeg(localSpikePhases))
                        disp(circ_r(ang2rad(localSpikePhases)))
                        disp('')
                
                    end
                    %}
                    if(circ_r(ang2rad(localSpikePhases(:)))>0) %consistency threshold, otherwise no mean phase defined
                        localAvgPhase=circMeanDeg(localSpikePhases);
                        %localAvgPhase=nanmean(localSpikePhases); %not so
                        %circular on a theta time scale.... but offset
                        %could be due to wrong source LFP
                        localAvgPhasePerTime(ti)=localAvgPhase;
                    %else
                    %    localSpikePhases
                    %    disp('')
                        
                    end
                    
                end
                
                 if(~isempty(localSpikeTimes))
                    localAvgTiming=nanmean(localSpikeTimes);
                    localAvgTimingPerTime(ti)=localAvgTiming;
                end
        end

end


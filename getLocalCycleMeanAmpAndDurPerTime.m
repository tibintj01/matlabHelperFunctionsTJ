function [localCycleAmpPerTime,localCycleDurPerTime,localCycleNumPerTime,localCycleStartTimePerTime] =getLocalCycleMeanAmpPerTime(timeBinCenters,cycleStartTimes,cycleAmps,cycleDurs)
%gets local running circular avg phase matching timebins in timeBinCenters based on cycle times

        numTimeBins=length(timeBinCenters);

        approxTimeStep=median(diff(timeBinCenters));

        localAvgCycleAmpPerTime=NaN(size(timeBinCenters));
        
	for ti=1:numTimeBins
		currTime=timeBinCenters(ti);
		[~,closestCycleStartID]=min(abs(cycleStartTimes-currTime));
        
        %if closest is next, revert to last - so each start gets whole
        %cycle after it
        if(cycleStartTimes(closestCycleStartID)>currTime && closestCycleStartID>1)
            closestCycleStartID=closestCycleStartID-1;
        end

		closestCycleAmp=cycleAmps(closestCycleStartID);
		closestCycleDur=cycleDurs(closestCycleStartID);
        closestCycleNum=closestCycleStartID;
        closestCycleStartTime=cycleStartTimes(closestCycleStartID);
                
        localCycleAmpPerTime(ti)=closestCycleAmp;
        localCycleDurPerTime(ti)=closestCycleDur;
        localCycleNumPerTime(ti)=closestCycleNum;
        localCycleStartTimePerTime(ti)=closestCycleStartTime;
    end




function [localFiringRatePerTime]=getLocalWindowSpikeRate(spikeTimes,timeBinCenters,localWindowWidth)

	numTimeBins=length(timeBinCenters);

	approxTimeStep=median(diff(timeBinCenters));

	localWindowHalfWidthIdx=round((localWindowWidth/approxTimeStep)/2);

	localFiringRatePerTime=NaN(size(timeBinCenters));
	for ti=1:numTimeBins
		localWindStartIdx=max(1,ti-localWindowHalfWidthIdx);
		localWindEndIdx=min(numTimeBins,ti+localWindowHalfWidthIdx);

		localWindStartTime=timeBinCenters(localWindStartIdx);
		localWindEndTime=timeBinCenters(localWindEndIdx);
		localWindWidth=localWindEndTime-localWindStartTime;
        
		localSpikeCount=sum(spikeTimes>=localWindStartTime & spikeTimes <localWindEndTime);
        
        
       

		localSpikeRate=localSpikeCount/localWindWidth;
		localFiringRatePerTime(ti)=localSpikeRate;
	end
	

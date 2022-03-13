function [rateTimeBinCenters,firingRateOverTime] = spikeTimesToFiringRate(spikeTimes,binWidth,minTime,maxTime)
	if(~exist('maxTime'))
		maxTime =  max(spikeTimes) +  1;   % what is the time of the last spike. How else could you have coded this line?
	end
	
	if(~exist('minTime'))
		minTime=0;
	end
	rateTimeBinEdges =  minTime:binWidth:maxTime;
	rateTimeBinCenters = (minTime+binWidth/2):binWidth:(maxTime-binWidth/2);
	%spikeRate =  histc(spikeTimes, rateTimeBinEdges);
	spikeRate = histcounts(spikeTimes, rateTimeBinEdges);
	%spikeRate =  spikeRate(1:end-1);
	firingRateOverTime =  spikeRate/binWidth;	

	

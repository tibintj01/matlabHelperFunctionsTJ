function [rateTimeBinCenters,firingRateOverTime] = getFiringRateOverTime(spikeTimes,rateWindow,maxTime)
	if(~exist('maxTime'))
		maxTime =  max(spikeTimes) +  1;   % what is the time of the last spike. How else could you have coded this line?
	end
	
	rateTimeBinEdges =  0:rateWindow:maxTime;
	rateTimeBinCenters = rateWindow/2:rateWindow:maxTime-rateWindow/2;
	spikeRate =  histc(spikeTimes, rateTimeBinEdges);
	spikeRate =  spikeRate(1:end-1);
	firingRateOverTime =  spikeRate/rateWindow;	

	

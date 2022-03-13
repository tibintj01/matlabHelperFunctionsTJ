function [binCenters,firingRate]=getDynamicFiringRate(spikeTimes,binWidth,minTime,maxTime)
	%Tibin John, 9/13/17
	%Description: Computes firing rate over time for a given cell
	%Input: vector of spike times, desired time width of bin, max time value
	%Output: firing rates at each bin center, and bin centers (time axis)

	if(~exist('maxTime'))
		maxTime=max(spikeTimes);
	end
	if(~exist('minTime'))
		minTime=0;
	end

	binEdges=minTime:binWidth:maxTime;

	binCenters=binEdges(2:end)-binWidth/2;

	firingRate=zeros(size(binCenters));
	for binNum=1:length(binCenters)
		currBinCenter=binCenters(binNum);
		spikeIdxes=find(spikeTimes>currBinCenter-binWidth/2 &...
				 spikeTimes < currBinCenter+binWidth/2);
		spikeCount=length(spikeIdxes);

		firingRate(binNum)=spikeCount/binWidth;
	end



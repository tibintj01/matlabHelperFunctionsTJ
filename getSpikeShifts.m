function [shiftsToBest,alignableIdxes] = getTroughAlignedSpikeIdxes(spikeTroughIdxes)
	bestMinIdx=mode(spikeTroughIdxes);

	shiftsToBest=-(spikeTroughIdxes-bestMinIdx);

	alignableIdxes=find(abs(shiftsToBest)<2);
	 

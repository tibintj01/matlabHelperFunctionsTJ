function [isolatedSpikeIdxes]=getIsolatedSpikeTimes(spikeTimes,isoTimeBefore,isoTimeAfter)

	spikeTimes=spikeTimes(:);

	isiBefore=[0;diff(spikeTimes)];
	isiAfter=[diff(spikeTimes); 0];

	isolatedSpikeIdxes=find(isiBefore>isoTimeBefore & isiAfter > isoTimeAfter);



	

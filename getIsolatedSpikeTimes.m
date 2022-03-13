function [isolatedSpikeTimes]=getIsolatedSpikeTimes(spikeTimes,isoTimeBefore,isoTimeAfter)

	spikeTimes=spikeTimes(:);

	isiBefore=[0;diff(spikeTimes)]
	isiAfter=[diff(spikeTimes); 0]

	isolatedSpikeTimes=spikeTimes(isiBefore>isoTimeBefore & isiAfter > isoTimeAfter);



	

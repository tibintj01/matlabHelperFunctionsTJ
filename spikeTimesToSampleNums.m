function [sampleNums]=spikeTimesToSampleNums(spikeTimes,Fs)
	sampleNums=round(spikeTimes*Fs);

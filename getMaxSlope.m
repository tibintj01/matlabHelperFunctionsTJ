function [maxSlope]=getMaxSlope(waveform,Fs)

	interpFact=10;
        [interpolatedWaveform,newFs]=getInterped(waveform,interpFact,Fs);

        [minVal, minIdx]=min(interpolatedWaveform);
        [maxVal,maxIdx]=max(interpolatedWaveform);

	instDeriv=diff(interpolatedWaveform(minIdx:maxIdx));
	maxSlope=max(instDeriv)*newFs;

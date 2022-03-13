function [avgWaveform,sdWaveform]=getAvgWaveformGivenCell(cellStruct)
	waveformsMatrix=cellStruct.waveForms;
	avgWaveform=nanmean(waveformsMatrix,1);
	sdWaveform=nanstd(waveformsMatrix,0,1);

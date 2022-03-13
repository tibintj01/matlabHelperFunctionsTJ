function [currWaveformTPwidth,currWaveformTPamp]=getWidthAmpOfWaveform(currWaveform)
	saveStr='individualThreshCrossUnsorted';
	cellPropSaveDir=sprintf('/nfs/turbo/lsa-ojahmed/tibin/processedHumanData/seizurePaper/blackDotSpikeProps')
	sessionIdx=1;
	ch=1;
	%[waveformProps]=getCellPropGivenWaveforms(currWaveform,0,ch,sessionIdx,saveStr,cellPropSaveDir,1);
	[waveformProps]=getSingleWaveformProp(currWaveform,0,ch,sessionIdx,saveStr,cellPropSaveDir,1);

	currWaveformTPwidth=waveformProps.spikeTPWidth;
	currWaveformTPamp=waveformProps.spikeTPAmp;


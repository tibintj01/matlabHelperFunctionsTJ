function [filteredWaveforms] = getFilteredWaveforms(waveforms,Fs)
	%why does this code propagate NaN's? - possibly Omar's filterLFP

	filteredWaveforms=zeros(size(waveforms));

	disp('filtering each waveform......')
	tic

	wavesT=waveforms';
	concatWavs=wavesT(:);

	%filteredConcat=getFilteredLFPforSpikeSorting(waveform,Fs);
	filteredConcat=getFilteredLFPforSpikeSorting(concatWavs,Fs);

	filteredWaveforms=reshape(filteredConcat,[size(waveforms,2), size(waveforms,1)])';
	%for row=1:size(waveforms,1)
	%	waveform=waveforms(row,:);
	%	filteredWaveforms(row,:)=getFilteredLFPforSpikeSorting(waveform,Fs);
	%end
	toc

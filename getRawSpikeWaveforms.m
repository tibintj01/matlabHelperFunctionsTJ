function [rawSpikeWaveforms,newCell] = getRawSpikeWaveforms(cell)
	%Description: given cell 
	%Input:
	%Output:

	spikeTimes=cell.spikeTimes;
	Fs=cell.Fs;
	chanNum=cell.chanNum;

	totalLength=48;
	samplesBeforeCrossing=9;
	samplesAfterCrossing=totalLength-samplesBeforeCrossing-1;

	rawSpikeWaveforms=NaN(length(spikeTimes),totalLength);

	disp('loading concatenated file.......')
	tic
	lfpDataFile=loadConcatFile(chanNum);
	rawLFP=lfpDataFile.concatenatedLFP;
	toc
	disp('done')

	tic
	disp('getting raw segment for each spike....')
	for spikeNum=1:length(spikeTimes)
		spikeTime=spikeTimes(spikeNum);
		startTime=spikeTime-samplesBeforeCrossing/Fs;
		endTime=spikeTime+samplesAfterCrossing/Fs;
		%rawSpikeWaveform=getRawSegment_MG49_sess3(startTime,endTime,chanNum);
		%rawSpikeWaveform=getRawSegmentFromConcat(startTime,endTime,chanNum);
		startSample=round(startTime*Fs);
	        endSample=round(endTime*Fs);

		rawSpikeWaveform=rawLFP(startSample:endSample);
		%rawSpikeWaveform=lfpDataFile.concatenatedLFP(1,startSample:endSample);
		rawSpikeZeroed=rawSpikeWaveform-median(rawSpikeWaveform);
		rawSpikeWaveforms(spikeNum,1:length(rawSpikeZeroed))=rawSpikeZeroed(:)';
	end
	toc
	disp('done')
	newCell=cell;
	newCell.rawSpikeWaveforms=rawSpikeWaveforms;

	rawSpikeWaveforms=rawSpikeWaveforms';


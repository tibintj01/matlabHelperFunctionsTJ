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
	padFront=0.007;
	padEnd=0.007;

	totalLength=totalLength+padFront*Fs+padEnd*Fs;

	rawSpikeWaveforms=NaN(length(spikeTimes),totalLength);

	disp('loading concatenated file.......')
	tic
	%lfpDataFile=loadConcatFile(chanNum);
	%rawLFP=lfpDataFile.concatenatedLFP;
	%[recordingDir,startSample,endSample]=getFileAndSampleNumbers(chanNum,timeStart,timeEnd);
	%parentDirName='C:/Users/tibin/spikeDynamicsAnalysisTibin/session3Mat/MG49';
	%fullPath=fullfile(parentDirName,recordingDir,sprintf('lfp%d_original.mat',chanNum))
	%rawLFPdata=load(fullPath);
	%rawLFP=rawLFPdata.lfp;
	sessNum=3;
	dataDir='C:/Users/tibin/spikeDynamicsAnalysisTibin/session3Mat/';
	rawLFP=getConcatenatedLFPforChan(sessNum,dataDir,chanNum);
	toc
	disp('done')
	maxIdx=length(rawLFP);

	tic
	disp('getting extended raw segment for each spike....')
	for spikeNum=1:length(spikeTimes)
		spikeTime=spikeTimes(spikeNum);
		startTime=spikeTime-(samplesBeforeCrossing/Fs)-padFront;
		endTime=spikeTime+samplesAfterCrossing/Fs+padEnd;
		%rawSpikeWaveform=getRawSegment_MG49_sess3(startTime,endTime,chanNum);
		%rawSpikeWaveform=getRawSegmentFromConcat(startTime,endTime,chanNum);
		startSample=max(1,round(startTime*Fs));
	        endSample=min(round(endTime*Fs),maxIdx);
		spikeIdxInRaw=round(spikeTime*Fs)-startSample+1;

		rawSpikeWaveform=rawLFP(startSample:endSample);
		%rawSpikeWaveform=lfpDataFile.concatenatedLFP(1,startSample:endSample);
		rawSpikeZeroed=rawSpikeWaveform-median(rawSpikeWaveform);
		%rawSpikeWaveforms(spikeNum,1:length(rawSpikeZeroed))=rawSpikeZeroed(:)';
		rawSpikeWaveforms(spikeNum,1:length(rawSpikeWaveform))=rawSpikeWaveform(:)';
	end
	toc
	disp('done')
	newCell=cell;
	newCell.rawSpikeWaveforms=rawSpikeWaveforms;
	newCell.padFrontTime=padFront;
	newCell.padEndTime=padEnd;
	newCell.spikeIdxInRaw=spikeIdxInRaw;

	rawSpikeWaveforms=rawSpikeWaveforms';


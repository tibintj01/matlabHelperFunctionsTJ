function [rawSpikeWaveforms] = getRawWaveformsOfCell(ns5dirName,descriptiveFilename,ch,spikeTimes,samplesBeforeCrossing,samplesAfterCrossing)
%function [rawSpikeWaveformsPlusMinus2msMat,rawSpikeWaveforms] = getRawWaveformsOfCell(ns5dirName,descriptiveFilename,ch,spikeTimes,samplesBeforeCrossing,samplesAfterCrossing)

	ns5FileName=fullfile(ns5dirName,sprintf([descriptiveFilename '_ch%d.ns5'],ch));
	rawDataStruct=openNSx(ns5FileName,'read','report');

	Fs=rawDataStruct.MetaTags.SamplingFreq;

	rawLFP=double(rawDataStruct.Data);

	rawSpikeWaveformsPlusMinus2msMat=zeros(length(spikeTimes),samplesBeforeCrossing+samplesAfterCrossing+1+Fs*4);
	
	tic
        disp('getting raw segment for each spike....')
        for spikeNum=1:length(spikeTimes)
                spikeTime=spikeTimes(spikeNum);
                startTime=spikeTime-samplesBeforeCrossing/Fs;
		%startTimeMinus=max(0,startTime-2);

                endTime=spikeTime+samplesAfterCrossing/Fs;
		%endTimePlus=min(endTime+2,length(rawLFP)/Fs);
                %rawSpikeWaveform=getRawSegment_MG49_sess3(startTime,endTime,chanNum);
                %rawSpikeWaveform=getRawSegmentFromConcat(startTime,endTime,chanNum);
          
	      startSample=round(startTime*Fs);
                endSample=round(endTime*Fs);

	      %startSampleMinus=round(startTimeMinus*Fs);
                %endSamplePlus=round(endTimePlus*Fs);

                rawSpikeWaveform=rawLFP(startSample:endSample);
		%rawSpikeWaveformsPlusMinus2ms=rawLFP(startSampleMinus:endSamplePlus);
                %rawSpikeWaveform=lfpDataFile.concatenatedLFP(1,startSample:endSample);

                rawSpikeZeroed=rawSpikeWaveform-median(rawSpikeWaveform);
                rawSpikeWaveforms(spikeNum,1:length(rawSpikeZeroed))=rawSpikeZeroed(:)';

		%rawSpikeWaveformsPlusMinus2msMat(spikeNum, 1:length(rawSpikeWaveformsPlusMinus2ms))=rawSpikeWaveformsPlusMinus2ms(:)';

        end
        toc

	rawSpikeWaveforms=rawSpikeWaveforms';


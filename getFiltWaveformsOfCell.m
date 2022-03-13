function [filteredWaveforms] = getFiltWaveformsOfCell(ns5dirName,descriptiveFilename,ch,spikeTimes,samplesBeforeCrossing,samplesAfterCrossing)
%function [filtSpikeWaveformsPlusMinus2msMat,filtSpikeWaveforms] = getRawWaveformsOfCell(ns5dirName,descriptiveFilename,ch,spikeTimes,samplesBeforeCrossing,samplesAfterCrossing)
	plotRawWaveforms=0;

	ns5FileName=fullfile(ns5dirName,sprintf([descriptiveFilename '_ch%d.ns5'],ch));
	rawDataStruct=openNSx(ns5FileName,'read','report');

	Fs=rawDataStruct.MetaTags.SamplingFreq;


	rawLFP=double(rawDataStruct.Data);
	fds


	lfpTimeSeries.lfp=rawLFP;
	lfpTimeSeries.Fs=Fs;

	%filtLFP=getFilteredLFPforSpikeSorting(rawLFP,Fs);

	%rawSpikeWaveformsPlusMinus2msMat=zeros(length(spikeTimes),samplesBeforeCrossing+samplesAfterCrossing+1+Fs*4);
	%test
	deadSpace=6000; %- takes a few hours to run...
	[rawSpikeWaveforms] = getWaveformsFromLFP(lfpTimeSeries,spikeTimes,deadSpace);
	rawSpikeWaveforms=rawSpikeWaveforms';
	%{
	%deadSpace=500; %16.67ms on either side (5 cycles of 300 Hz to avoid edge effects) - enough?
	realFirstIdx=deadSpace+1;
	realLastIdx=realFirstIdx+samplesBeforeCrossing+1+samplesAfterCrossing-1;

	rawSpikeWaveforms=zeros(length(spikeTimes),deadSpace+48+deadSpace);
	%rawSpikeWaveforms=NaN(length(spikeTimes),deadSpace+48+deadSpace);
	
	tic
        disp('getting filt segment for each spike....')
	numWithinTimeSpikes=0;
	for spikeNum=1:length(spikeTimes)
                spikeTime=spikeTimes(spikeNum);
                startTime=spikeTime-samplesBeforeCrossing/Fs;
		%startTimeMinus=max(0,startTime-2);

                endTime=spikeTime+samplesAfterCrossing/Fs;
		%endTimePlus=min(endTime+2,length(filtLFP)/Fs);
                %filtSpikeWaveform=getRawSegment_MG49_sess3(startTime,endTime,chanNum);
                %filtSpikeWaveform=getRawSegmentFromConcat(startTime,endTime,chanNum);
          
	        startSample=max(1,round(startTime*Fs)-deadSpace);
                endSample=min(length(rawLFP),round(endTime*Fs)+deadSpace);

	        %startSampleMinus=round(startTimeMinus*Fs);
                %endSamplePlus=round(endTimePlus*F

		if(startSample > length(rawLFP) || endSample > length(rawLFP))
			continue
		end
		%fds
		numWithinTimeSpikes=numWithinTimeSpikes+1;

                rawSpikeWaveform=rawLFP(startSample:endSample);
		rawSpikeZeroed=rawSpikeWaveform-median(rawSpikeWaveform(realFirstIdx:realLastIdx));
               %rawSpikeZeroed=rawSpikeWaveform;

		 rawSpikeWaveforms(spikeNum,1:length(rawSpikeWaveform))=rawSpikeZeroed(:)';

		%filtSpikeWaveformsPlusMinus2msMat(spikeNum, 1:length(filtSpikeWaveformsPlusMinus2ms))=filtSpikeWaveformsPlusMinus2ms(:)';

        end
        toc
	%}

	if(plotRawWaveforms)
                figure
                subplot(4,1,1)
                spikeNumsDisp=find(spikeTimes>100 & spikeTimes<500);
                plot(rawSpikeWaveforms(spikeNumsDisp,deadSpace:(deadSpace+48))')
                title('spikeTimes>100 & spikeTimes<300')
		ylim([-250 250])

                subplot(4,1,2)
                spikeNumsDisp=find(spikeTimes>500 & spikeTimes<900);
                plot(rawSpikeWaveforms(spikeNumsDisp,deadSpace:(deadSpace+48))')
                title('spikeTimes>500 & spikeTimes<700')
		ylim([-250 250])

                subplot(4,1,3)
                spikeNumsDisp=find(spikeTimes>1000 & spikeTimes<1400);
                plot(rawSpikeWaveforms(spikeNumsDisp,deadSpace:(deadSpace+48))')
                title('spikeTimes>1000 & spikeTimes<1200')
		ylim([-250 250])

                subplot(4,1,4)
                spikeNumsDisp=find(spikeTimes>2500 & spikeTimes<2900);
                plot(rawSpikeWaveforms(spikeNumsDisp,deadSpace:(deadSpace+48))')
                title('spikeTimes>2500 & spikeTimes<2700')
		ylim([-250 250])
	
                %uberTitle(sprintf('Human nex spike times in raw ns5 %s Ch %da','MG63',ch))
                maxFig
                saveas(gcf,'RawSpikeWaveformsTestTimesHumanData.tif')
                fds
        end

	[filteredWaveforms] = getFilteredWaveforms(rawSpikeWaveforms,Fs);
	filteredWaveforms=filteredWaveforms';
	filteredWaveforms=filteredWaveforms(realFirstIdx:realLastIdx,:);
	filteredWaveforms=filteredWaveforms(:,1:numWithinTimeSpikes);


	%fds
	%figure
	%subplot(1,2,1)
	%plot(filteredWaveforms(:,1:10))  
	%subplot(1,2,2)
	%rawDisplay=rawSpikeWaveforms(1:10,realFirstIdx:realLastIdx)';
	%plot(rawDisplay) 
	%ylim([-1000 600])
	%saveas(gcf,'testFiltWaves.tif') 
	%fds

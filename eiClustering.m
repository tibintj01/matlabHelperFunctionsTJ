function [] = eiClustering(cells,useRaw)

	figure

	for cellIdx=1:length(cells)
		disp(sprintf('%.4f complete.....',cellIdx/length(cells)))
		cell=cells{cellIdx};
		
		if(cell.qualityRating<3)
			continue
		end

		Fs=cell.Fs;

		%cellAvgWave=getAvgWaveform(cell);
		if(useRaw)
			waveforms=cell.rawSpikeWaveforms;
			startIdx=cell.spikeIdxInRaw-9;
			endIdx=cell.spikeIdxInRaw+38;
			%filtWaveforms=getFilteredWaveforms(waveforms(:,startIdx:endIdx),Fs);
			filtWaveforms=getFilteredWaveforms(waveforms,Fs);
			waveforms=filtWaveforms(:,startIdx:endIdx);
			%waveforms=filtWaveforms;
		
			%waveforms=getFilteredLFPforSpikeSorting(waveforms,Fs);
			%[waveforms]=getRawSpikeWaveforms(cell);
		else
			%online waveforms
			waveforms=cell.waveForms;
		end

		
		[meanVToP,meanHalfTrough]=getMeanShapeParameters(waveforms,Fs);	
		%[meanVToP,meanHalfTrough]=getParametersOnMeanWave(waveforms,Fs);	

		%avgPeriod=1/getAvgFiringRate(cell);
		plot(meanVToP*1000,meanHalfTrough*1000,'ko','MarkerSize',1)
		%plot(maxSlope,halfPeak*1000,'ko','MarkerSize',1)

		%plot3(avgPeriod*1000,vToP*1000,halfPeak*1000,'ko','MarkerSize',1)
		daspect([1 1 1])
		hold on
	end
	
	xlabel('Valley to Peak (ms)')
	%xlabel('Max slope (mV/sec)')
	ylabel('Half peak width (ms)')
	%zlabel('Avg ISI (ms)')
	%ylabel('Valley to Peak (ms)')
	%zlabel('Half peak width (ms)')
	%xlabel('Avg ISI (ms)')
	%xlim([0 10000])

function [] = plotCellsInWaveShapeSpace(cellWaveformMat)

	figure

	for cellIdx=1:length(cellWaveformMat)
		cell=cellWaveformMat{cellIdx};
		
		if(cell.qualityRating<3)
			continue
		end

		Fs=cell.Fs;

		cellAvgWave=getAvgWaveform(cell);
		vToP=getValleyToPeakTime(cellAvgWave,Fs);
		maxSlope=getMaxSlope(cellAvgWave,Fs);
		halfPeak=getHalfPeakWidth(cellAvgWave,Fs);
		
		avgPeriod=1/getAvgFiringRate(cell);
		plot(vToP*1000,halfPeak*1000,'ko','MarkerSize',1)
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

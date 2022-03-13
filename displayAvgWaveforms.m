function [] = displayAvgWaveforms(cellList)

	figure
	for cellIdx=1:length(cellList)
		cell=cellList{cellIdx};
		Fs=cell.Fs;
		avgWave=getAvgWaveform(cell);
		timeAxis=(1:length(avgWave))/Fs;

		plot(timeAxis,avgWave,'r-')
		hold on
		xlabel('Time (ms)')	
	end

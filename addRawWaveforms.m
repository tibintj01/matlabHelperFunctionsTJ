function [newCells]=addRawWaveforms(cells,cellSaveName)
	disp('adding raw waveforms.....')
	newCells=cell(size(cells));

	for cellIdx=1:length(cells)
		cellIdx/length(cells)
		currCell=cells{cellIdx};
		%[dummyVar,newCell]=getRawSpikeWaveforms(currCell);
		[dummyVar,newCell]=getExtendedRawSpikeWaveforms(currCell);
		%newCell
		newCells{cellIdx}=newCell;
	end

	disp('saving cells with raw waveforms added.....')
	tic
	save(cellSaveName,'newCells','-v7.3')
	toc
	disp('done')


function [allCellAvgWaveform]=getAvgWaveformAcrossCells(cellDir,szID)

	%szID='MG49-seizure45';

	filePaths=getRegexFilePaths(cellDir,sprintf('*cell_properties*%s*.mat',szID));
	
	numCells=length(filePaths);
	numBins=471; %result of splineDt=0.1 on 48 bin original waveform
	sumWaveforms=zeros(numBins,1);

	numCellsUsed=0;
	numSpikesUsed=0;
	for i=1:length(filePaths)
		i/numCells
		%avgWaveform=getAvgGoodWaveformForCell(filePaths{i});
		%sumWaveforms=nansum(cat(3,sumWaveforms(:),avgWaveform(:)),3);
		cellWaveforms=getGoodWaveforms(filePaths{i});
		if(length(cellWaveforms)>1)
			for c=1:size(cellWaveforms,2)
				sumWaveforms=nansum(cat(3,sumWaveforms(:),cellWaveforms(:,c)),3);
				%numCellsUsed=numCellsUsed+1;
				numSpikesUsed=numSpikesUsed+1;
			end
		end
	end

	%allCellAvgWaveform=sumWaveforms/numCells;
	allCellAvgWaveform=sumWaveforms/numSpikesUsed;

	save(fullfile(cellDir,sprintf('allCellAvgWaveform%s.mat',szID)),'allCellAvgWaveform')

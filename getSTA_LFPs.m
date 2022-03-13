function [] = getSTA_LFPs(cellPropDir)
	%inputCellArray=cell();

	freqBands=freqBandNames();
	cellPropFileNames=getRegexFilePaths(cellPropDir,'*cell_prop*mat');
	maxNumSpikes=1000;
	count=0;
	for i=1:length(freqBands)
		for j=1:length(cellPropFileNames)
			%count=count+1;
			%inputCellArray{count,1}=cellPropFileNames{j};
			%inputCellArray{count,2}=freqBands{i};
			%inputCellArray{count,3}=maxNumSpikes;
			plotPeriSpikeForCell_flux(cellPropFileNames{j},freqBands{i},maxNumSpikes)
		end
	end	
	%inputCellArray
	%save('LFP_sta_ParInput.mat','inputCellArray')

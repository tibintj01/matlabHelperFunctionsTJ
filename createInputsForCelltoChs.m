function [] = createInputsForCelltoChs(cellPropDir)
	filePaths=getRegexFilePaths(cellPropDir,'*cell_propert*mat');

	freqBands={'Delta','Alpha','LowGamma','HighGamma'};
	%freqBands={'Delta'}

	inputCellArray=cell(1,3);

	maxNumSpikes=10000000;
	inputCount=1;
	%for i=1:length(filePaths)
	for i=1:40
		for j=1:length(freqBands)
			inputCellArray{inputCount,1}=filePaths{i};
			inputCellArray{inputCount,2}=freqBands{j};
			inputCellArray{inputCount,3}=maxNumSpikes;

			inputCount=inputCount+1;	
		end
	end 

	inputCellArray
	save('inputCellToAllChannelsPhase.mat','inputCellArray')

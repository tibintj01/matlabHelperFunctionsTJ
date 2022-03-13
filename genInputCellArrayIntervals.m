function [] = genInputCellArrayIntervals(maxIdx,stepSize)
	inputCellArray=cell(1,2);

	count=1;
	for i=1:stepSize:maxIdx
		inputCellArray{count,1}=i;
		inputCellArray{count,2}=i+stepSize-1;
		count=count+1;
	end

	save('IdxIntervalsInputCell.mat','inputCellArray')

	 

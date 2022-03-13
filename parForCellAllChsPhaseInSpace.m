
inputPath='/nfs/turbo/lsa-ojahmed/inputCellToAllChannelsPhase.mat'

%poolobj=gcp
%addAttachedFiles(poolobj,inputPath)
%poolobj

data=load(inputPath);
numRuns=size(data.inputCellArray,1)

%parpool(20)
parpool(16)
parfor i= 1:numRuns
	data=load(inputPath);
	i
	data.inputCellArray{i,:}
	tic
	try
	getCellToAllChsSpatialPhase(data.inputCellArray{i,:});
	%getCellToAllChsSpatialPhaseSubset(data.inputCellArray{i,:});
	catch
		disp('skipping in parfor....')
	end
	toc
end

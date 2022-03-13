function [] = parForate(funtionPath,inputPath)

%poolobj=gcp
%addAttachedFiles(poolobj,fileName)
%poolobj
data=load(inputPath);
numRuns=size(data.inputCellArray,1);

parpool(20)
parfor i= 1:numRuns
	data=load(inputPath);
	tic
	decimateConcatenedNS5s(data.inputCellArray{i,:});
	toc
end

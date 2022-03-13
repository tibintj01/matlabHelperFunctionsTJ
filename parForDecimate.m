
inputPath='/nfs/turbo/lsa-ojahmed/mg29decimateInputs.mat';

%poolobj=gcp
%addAttachedFiles(poolobj,inputPath)
%poolobj

data=load(inputPath);
numRuns=size(data.inputCellArray,1)

try
        parpool(20)
catch
        display('new worker pool not activated....')
end
tic
%parfor i= 1:numRuns
parfor i= 1:20
	data=load(inputPath);
	%tic
	decimateConcatenedNS5s(data.inputCellArray{i,:});
	%toc
end
toc

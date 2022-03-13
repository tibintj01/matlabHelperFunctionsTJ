function [] = runMatlabParInst(functionName,inputMatPath,parIdx)

	funcHandle=str2func(functionName);

	inputMat=load(inputMatPath);

	fieldNames=fieldnames(inputMat);

	inputCellArray=inputMat.(fieldNames{1});

	instanceInput=inputCellArray(parIdx,:);

	funcHandle(instanceInput{:})
			

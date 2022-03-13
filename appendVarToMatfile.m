function [] = appendVarToMatfile(matFilePath,value,varName)
	
	save(matFilePath,varName,

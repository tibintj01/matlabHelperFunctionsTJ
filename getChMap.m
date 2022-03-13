function [arrayMap electrodeXY] = getChMap(ptDir)
	[arrayMap electrodeXY electrodeImp] = neuroportArrayData(ptDir);
	arrayMap=fliplr(arrayMap)';


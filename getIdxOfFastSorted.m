function [I] = getIdxOfFastSorted(desiredValues,values)
	edges = [-Inf, mean([values(2:end); values(1:end-1)]), +Inf];
	I = discretize(desiredValues, edges);

function [szIdx] = getSzIdx(propPath)
	%saved as seizure%d.mat
	szIdx=str2num(propPath(end-4));

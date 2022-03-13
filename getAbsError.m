function [absError] = getAbsError(vals,truth)
	absError=sum(abs(vals-truth));

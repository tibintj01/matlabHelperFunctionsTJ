function [zScored]=nanzscore(vals)
	zScored=(vals-nanmean(vals))/nanstd(vals);

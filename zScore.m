function [zScored]=zScore(vals)
	zScored=(vals-nanmean(vals))/nanstd(vals);

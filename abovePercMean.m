function [topMean] = topPercMean(values,topPerc)
	topFloor=prctile(values,topPerc);

	topMean=nanmean(values(values>topFloor));


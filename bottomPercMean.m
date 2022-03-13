function [bottomMean] = bottomPercMean(values,bottomPerc)
	bottomCeiling=prctile(values,bottomPerc);

	bottomMean=nanmean(values(values<bottomCeiling));


function [avgFiringRate]=getAvgFiringRate(cell)
	totalTime=cell.tend-cell.tbeg;

	totalNumSpikes=length(cell.spikeTimes);

	avgFiringRate=totalNumSpikes/totalTime;

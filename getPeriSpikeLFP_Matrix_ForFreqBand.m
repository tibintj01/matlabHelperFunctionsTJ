function [periTime,periSpikeLFP,periSEM,numSpikesUsed] = getPeriSpikeLFP_Matrix_ForFreqBand(spikeTimes,lfp,lfpTime,lowFreq,highFreq)

	lfpFs=round(1/median(diff(lfpTime)));
	%filter order, zscore
	filteredLFP=filterLFP(lfp,lfpFs,lowFreq,highFreq,20,1);

	preTime=1/lowFreq;
	postTime=1/lowFreq;

	spikeTimes=spikeTimes(spikeTimes>min(lfpTime) & spikeTimes <max(lfpTime));

	disp('starting getIdxOf.....')
	tic

	%spikeBins=getIdxOf(spikeTimes,lfpTime);
	%spikeBins=round(spikeTimes*lfpFs);
	spikeBins=getIdxOfFastSorted(spikeTimes,lfpTime);
	toc
	disp('done')
	normAmpAroundSpike=0;

	disp('starting omar function')
	tic
	%[periSpikeLFP periTime periSEM numSpikesUsed spikeIdxesUsed]=omarGetAvgLFP_Tibin(filteredLFP,lfpFs,spikeBins,-preTime,postTime,normAmpAroundSpike);
	[periSpikeLFP periTime periSEM numSpikesUsed spikeIdxesUsed]=omarGetPeriLFPMatrix_Tibin(filteredLFP,lfpFs,spikeBins,-preTime,postTime,normAmpAroundSpike);
	toc


function [avgLFP] = lfpAvgAnon(spikeIdxes,lfp)

	disp('in anon function')
	preIdx=30000/15;
	postIdx=30000/15;
	spikeIdxes=spikeIdxes(spikeIdxes>preIdx & spikeIdxes<length(lfp)-postIdx);
	getPeriLFPs=@(x) lfp((x-preIdx):(x+postIdx));

	periLFPcell=arrayfun(getPeriLFPs,spikeIdxes,'UniformOutput',false)
	periLFPs=cell2mat(periLFPcell(:));

	avgLFP=mean(periLFPs,1)

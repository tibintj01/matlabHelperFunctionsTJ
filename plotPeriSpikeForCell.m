function [] = plotPeriSpikeForCell(cell,titleInfo,maxNumSpike)

	[lfp,lfpFs]=getDecLFPforCell(cell);

	lfpTime=(1:length(lfp))/lfpFs;
	spikeTimes=cell.spikeTimes;

	spikeIdxes=1:(min(maxNumSpike,length(spikeTimes)));
	spikeTimes=spikeTimes(spikeIdxes);

	lowFreq=5;
	
	highFreq=12;
	
	disp('getting spike triggered average.......')
	tic
	[periTime,periSpikeLFP,periSEM,numEvents]=getPeriSpikeLFPavg(spikeTimes,lfp,lfpTime,lowFreq,highFreq);
	%[avgLFP]= testAnonymous(getIdxOf(spikeTimes,lfpTime),lfp);
	%[avgLFP]= testAnonymous(randi(length(lfp),size(spikeTimes)),lfp);
	toc

	figure
	if(lowFreq>=5)
	shadedErrorBar(periTime*1000,periSpikeLFP,periSEM)
	xlabel('Time after spike (msec)')
	else
	shadedErrorBar(periTime*1,periSpikeLFP,periSEM)
	xlabel('Time after spike (sec)')
		
	end

	ylabel('LFP')

	title(sprintf([titleInfo ', Avg LFP (n=%d) Filtered from %d to %d Hz'],numEvents,lowFreq,highFreq))
	saveas(gcf,['avgLFP-' titleInfo '.tif'])

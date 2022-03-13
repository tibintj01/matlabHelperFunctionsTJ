function [] = saveAllLFPsingleCyclesForFreqBand(ch,freqBand,saveDir)
	if(exist('maxNumSpike','var') &&  isstr(maxNumSpike))
		maxNumSpike=str2num(maxNumSpike);
	end


	%get decimated LFP from this cell's channel
	[lfp,lfpFs]=getDecLFPforCh(ch);

	lfpTime=(1:length(lfp))/lfpFs;

	%get frequency boundaries corresponding to this band (e.g. 'Beta')
	[lowFreq,highFreq]=getFreqBand(freqBand);
	[lowFastFreq,highFastFreq]=getFastFreqBand(freqBand);
	%lowFreq=5;
	%highFreq=12;
	
	%compute average peri-spike LFP and/or phase statistics for this cell 
	disp('getting single cycles.......')
	tic
	%[periTime,periSpikeLFP,periSEM,numEvents]=getPeriSpikeLFPavg(spikeTimes,lfp,lfpTime,lowFreq,highFreq);
	%[lfpPhaseHist]=getLFPspikePhaseHist(spikeTimes,lfp,lfpTime,lowFreq,highFreq);
	chStr=num2str(ch)
	if(ch<10)
		chStr=['0' chStr];
	end
	%[singleCycleInfo]=humanProcessLFPMinMax(lfpTime,lfp,lfpFs,lowFreq,highFreq,['Ch' chStr '-singleCycle-' freqBand],1,0,1,1,1,1); 
	%[singleCycleInfo]=humanProcessLFPMinMax(lfpTime,lfp,lfpFs,lowFreq,highFreq,['Ch' chStr '-singleCycle-' freqBand],1,0,1,1,1,1); 
	centerOnMin=1;
	asymData = omarAsymmetryLFP(lfpTime, lfp, lfpFs, lowFreq, highFreq, lowFastFreq, highFastFreq, ['Ch' chStr '-singleCycle-' freqBand], centerOnMin);
	toc
	
	chStr=num2str(ch);

	if(ch<10)
		chStr=['0' chStr];
	end

	save(fullfile(saveDir,[freqBand 'SingleCycleBroadbandStats-Ch' chStr '.mat']),'ch','asymData','lowFreq','highFreq','lowFastFreq','highFastFreq')




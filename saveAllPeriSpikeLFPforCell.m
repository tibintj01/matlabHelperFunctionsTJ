function [] = saveAllPeriSpikeLFPforCell(cellPropFile,freqBand,maxNumSpike,saveDir)
	if(exist('maxNumSpike','var') &&  isstr(maxNumSpike))
		maxNumSpike=str2num(maxNumSpike);
	end

	cellProp=load(cellPropFile);

	%get phenotype of this cell
	if(cellProp.isInterneuronCell==0)
		cellType='RS';
	elseif(cellProp.isInterneuronCell==1)
		cellType='FS';
	else
		cellType='GS';
	end

	%get decimated LFP from this cell's channel
	[lfp,lfpFs]=getDecLFPforCell_flux(cellPropFile);

	lfpTime=(1:length(lfp))/lfpFs;

	%use only good (single, trough at expected time relative to threshold) spikes
	spikeTimes=cellProp.spikeTimes(cellProp.goodSpikes);


	%get only first N spikes if max is specified for quick run
	if(exist('maxNumSpike','var'))
		spikeIdxes=1:(min(maxNumSpike,length(spikeTimes)));
		spikeTimes=spikeTimes(spikeIdxes);
	end

	%get frequency boundaries corresponding to this band (e.g. 'Beta')
	[lowFreq,highFreq]=getFreqBand(freqBand);
	%lowFreq=5;
	%highFreq=12;
	
	%compute average peri-spike LFP and/or phase statistics for this cell 
	disp('getting spike triggered average.......')
	tic
	%[periTime,periSpikeLFP,periSEM,numEvents]=getPeriSpikeLFPavg(spikeTimes,lfp,lfpTime,lowFreq,highFreq);
	%[lfpPhaseHist]=getLFPspikePhaseHist(spikeTimes,lfp,lfpTime,lowFreq,highFreq);
	[periTime,periSpikeLFPMatrix,periSEM,numEvents]=getPeriSpikeLFP_Matrix_ForFreqBand(spikeTimes,lfp,lfpTime,lowFreq,highFreq);
	
	toc


	titleInfo=getFileNameFromPath(cellPropFile);
                titleInfo=titleInfo(1:(end-4));
                titleInfo=strrep(titleInfo,'_','-');
                titleInfo=strrep(titleInfo,'-cell-properties','');
                titleString=sprintf([cellType '-' freqBand '-' titleInfo ', LFP (n=%d) Filtered from %d to %d Hz'],numEvents,lowFreq,highFreq)

	save(fullfile(saveDir,[freqBand 'PeriSpikeLFPs-' titleInfo '-' cellType '.mat']),'periSpikeLFPMatrix','periTime','titleString','lowFreq','highFreq')



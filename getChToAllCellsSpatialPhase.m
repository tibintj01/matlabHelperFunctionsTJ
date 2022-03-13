function [] = getLFP_SpikePhaseStatsCellProp(cellPropFile,freqBand,maxNumSpike)
	if(exist('maxNumSpike','var') &&  isstr(maxNumSpike))
		maxNumSpike=str2num(maxNumSpike);
	end

	cellPropDir=getDirNameFromPath(cellPropFile);
	allCellPropPaths=getRegexFilePaths(cellPropDir,'*cell_propert*mat');	

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


	%get frequency boundaries corresponding to this band (e.g. 'Beta')
	[lowFreq,highFreq]=getFreqBand(freqBand);
	%lowFreq=5;
	%highFreq=12;
	
		titleInfo=getFileNameFromPath(cellPropFile);
		titleInfo=titleInfo(1:(end-4));
		titleInfo=strrep(titleInfo,'_','-');
		titleInfo=strrep(titleInfo,'-cell-properties','');
	%compute average peri-spike LFP and/or phase statistics for this cell 
	disp('getting spike-lfp phase stats.......')
	tic
	phaseHistTitleStr=[titleInfo '-' cellType];
	% '-' freqBand];
	%[periTime,periSpikeLFP,periSEM,numEvents]=getPeriSpikeLFPavg(spikeTimes,lfp,lfpTime,lowFreq,highFreq);
	%[lfpPhaseHist]=getLFPspikePhaseHist(spikeTimes,lfp,lfpTime,lowFreq,highFreq);
	close all
	figure(10);
	%clrmp=copper(80);
	%colormap(copper(clrmp))
       	%colorbar('Units','normalized','Position',  
	%ha = tight_subplot(10,10,[.01 .03],[.1 .01],[.01 .01]);
        ha = tight_subplot(10,10,[.01 .03],[.1 .1],[.01 .1]);
                %axes('Color','none','XColor','none','YColor','none')
	%axes(ha(1))
       	%ax1Pos=get(ha(10),'Position')
	colormap(copper)
	%position is from left, from bottom, width,
	%cHandle=colorbar('Position',[0.95 0.03 0.007 0.9]) 
	cHandle=colorbar('Position',[0.91 0.1 0.007 0.8]) 
	set(cHandle,'FontSize',8) 
	maxKappa=0.7;
	%maxKappa=2;
	caxis([0 maxKappa])
	ylabel(cHandle,'\kappa','FontSize',16)
	%saveas(gcf,'testColorbar.tif')
	%fds
	 suptitle(['RS Phase Histograms, locking to LFP at ' phaseHistTitleStr ', ' freqBand])
          %set(get(gca,'title'),'Units','normalized','Position',[0.5 -0.1])
	
	%%loop through cell file names
	for cellPropPathNum=1:length(allCellPropPaths)
		cellProp=load(allCellPropPaths{cellPropPathNum});
	
		if(cellProp.isInterneuronCell ~=0 )
			continue
		end
		%use only good (single, trough at expected time relative to threshold) spikes
		spikeTimes=cellProp.spikeTimes(cellProp.goodSpikes);
		if(length(spikeTimes)<50)
			continue
		end
		%get only first N spikes if max is specified for quick run
		if(exist('maxNumSpike','var'))
			spikeIdxes=1:(min(maxNumSpike,length(spikeTimes)));
			spikeTimes=spikeTimes(spikeIdxes);
		end
		allCellPropPaths{cellPropPathNum}	
		[row,col]=getSpatialRowColFromFile(allCellPropPaths{cellPropPathNum})

		%figPos=sub2ind([10 10],row,col)
		figPos=sub2ind([10 10],col,row)
		[spikePhaseMean,spikePhaseKappa]=getLFP_SpikePhaseStats(spikeTimes,lfp,lfpTime,lowFreq,highFreq,[phaseHistTitleStr '-' freqBand],ha,figPos,maxKappa);
		%set(ha(figPos),'FontSize',3)

	end

	toc

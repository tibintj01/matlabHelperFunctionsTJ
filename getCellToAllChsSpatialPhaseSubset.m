function [] = getCellToAllChsSpatialPhaseSubset(cellPropFile,freqBand,maxNumSpike)
	if(exist('maxNumSpike','var') &&  isstr(maxNumSpike))
		maxNumSpike=str2num(maxNumSpike);
	end

	cellPropDir=getDirNameFromPath(cellPropFile);
	allCellPropPaths=getRegexFilePaths(cellPropDir,'*cell_propert*mat');	

	[cellRow,cellCol]=getSpatialRowColFromFile(cellPropFile)
	cellProp=load(cellPropFile);

	%get phenotype of this cell
	if(cellProp.isInterneuronCell==0)
		cellType='RS';
	elseif(cellProp.isInterneuronCell==1)
		cellType='FS';
	else
		cellType='GS';
	end



	%get frequency boundaries corresponding to this band (e.g. 'Beta')
	[lowFreq,highFreq]=getFreqBand(freqBand);
	%lowFreq=5;
	%highFreq=12;
	cellPropFileName=getFileNameFromPath(cellPropFile);
	
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
        %ha = tight_subplot(10,10,[.01 .03],[.1 .1],[.01 .1]);
	electrodeDiameter=5;
        %ha = tight_subplot(electrodeDiameter,electrodeDiameter,[.01 .03],[.1 .1],[.01 .1]);
        ha = tight_subplot(electrodeDiameter,electrodeDiameter,[.01 .01],[.1 .1],[.01 .1]);
                %axes('Color','none','XColor','none','YColor','none')
	%axes(ha(1))
       	%ax1Pos=get(ha(10),'Position')
	colormap(copper)
	%position is from left, from bottom, width,
	%cHandle=colorbar('Position',[0.95 0.03 0.007 0.9]) 
	cHandle=colorbar('Position',[0.91 0.1 0.007 0.8]) 
	set(cHandle,'FontSize',8) 
	if(contains(freqBand,'Delta'))
		maxKappa=2;
	else
		maxKappa=0.7;
	end

	caxis([0 maxKappa])
	ylabel(cHandle,'\kappa','FontSize',16)
	%saveas(gcf,'testColorbar.tif')
	%fds
	 suptitle(['RS Phase Histograms, locking to LFP at ' phaseHistTitleStr ', ' freqBand])
          %set(get(gca,'title'),'Units','normalized','Position',[0.5 -0.1])

	
	spikeTimes=cellProp.spikeTimes(cellProp.goodSpikes);
	lfpAlreadyDone=zeros(10,10);	
	cellSpatialPhaseLocking=zeros(10,10);
	

	localChNum=getChFromPath(cellPropFile);
	localChanNums=getLocalChNums(cellPropFileName,electrodeDiameter)

		
	%%loop through cell file names
	%for cellPropPathNum=1:length(allCellPropPaths)
	%for chanNum=1:96
	localColIdx=1;
	localRowIdx=1;
	localCount=1;
	%for chanNum=localChanNums
	for localRowIdx=1:electrodeDiameter
		for localColIdx=1:electrodeDiameter
			chanNum=localChanNums(localCount);
			localCount=localCount+1;
		%[row,col]=getSpatialRowColFromFile(allCellPropPaths{cellPropPathNum});
		%try
			try
			[row,col]=getSpatialRowColFromCh(cellPropFile,chanNum);
			catch
				disp('data not available....')
				continue
			end
			disp('getting lfp-cell relationship....')
			tic
			if(lfpAlreadyDone(row,col)==1)
				disp('this channel already calculated, continuing...')
				continue
			end

		
			%get decimated LFP from this cell's channel
			%[lfp,lfpFs]=getDecLFPforCell_flux(allCellPropPaths{cellPropPathNum});
			[lfp,lfpFs]=getDecLFPforCh(chanNum);
			lfpTime=(1:length(lfp))/lfpFs;
		
			%get only first N spikes if max is specified for quick run
			if(exist('maxNumSpike','var'))
				spikeIdxes=1:(min(maxNumSpike,length(spikeTimes)));
				spikeTimes=spikeTimes(spikeIdxes);
			end

			%figPos=sub2ind([10 10],row,col)
			%figPos=sub2ind([10 10],col,row)
			%figPos=sub2ind([electrodeDiameter electrodeDiameter],col,row)
			figPos=sub2ind([electrodeDiameter electrodeDiameter],localColIdx,localRowIdx)
			%[spikePhaseMean,spikePhaseKappa,spikePhases]=getLFP_SpikePhaseStats(spikeTimes,lfp,lfpTime,lowFreq,highFreq,[phaseHistTitleStr '-' freqBand],ha,figPos,maxKappa);
			[spikePhaseMean,spikePhaseKappa]=getLFP_SpikePhaseStats(spikeTimes,lfp,lfpTime,lowFreq,highFreq,[phaseHistTitleStr '-' freqBand],ha,figPos,maxKappa);
			%set(ha(figPos),'FontSize',3)
			lfpAlreadyDone(row,col)=1;
			
			cellSpatialMeanPhase(row,col)=spikePhaseMean;
			cellSpatialPhaseLocking(row,col)=spikePhaseKappa;
			toc
		%catch
		%	disp(sprintf('error in ch %d, skipping...',chanNum))
		%end
		end
	end
	%saveas(gcf,['PhaseHistogram-' titleStr '.tif'])
	%saveas(gcf,['PhaseHistogram-' titleStr '.tif'])
	print(gcf,['PhaseHistogram-' titleStr],'-dtiff','-r600')
	print(gcf,['PhaseHistogram-' titleStr],'-depsc','-r600')
	print(gcf,['PhaseHistogram-' titleStr],'-dpdf','-r600')
	%save(sprintf('%s-%s-SpatialPhaseAndLocking.mat',phaseHistTitleStr,freqBand),'cellSpatialMeanPhase','cellSpatialPhaseLocking','cellType','freqBand','cellRow','cellCol')

	toc

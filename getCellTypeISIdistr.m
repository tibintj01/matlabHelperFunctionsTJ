function [] = classifiedCellPropsDynamics(cellPropDir,groupName)
	%make sure input directory mat files are all same group
	%cellPropFilePaths=getFilePathsRegex(cellPropDir,'*_seizure1_*.mat')
	cellPropFilePaths=getFilePathsRegex(cellPropDir,'*.mat')
	%cellPropFilePaths=getFilePathsRegex(cellPropDir,'*MG49*_seizure1_*.mat')
	%cellPropFilePaths=getFilePathsRegex(cellPropDir,'*.mat')

	rsCount=0;
	fsCount=0;
	grayCount=0;
	 isi_RS=[];	
	 isi_FS=[];
	 isi_Gray=[];
	
	tic
	disp('getting firing rates by cell type.....')
	%save firing rate distributions by cell type	
	for cellIdx=1:length(cellPropFilePaths)
		cellProp=load(cellPropFilePaths{cellIdx});
		ISIvec=diff(cellProp.spikeTimes(cellProp.goodSpikes));
		
		if(cellProp.isInterneuronCell==0)
			rsCount=rsCount+1;
			%isi_RS(rsCount)=nanmean(1./diff(cellProp.spikeTimes));
			isi_RS=[isi_RS ISIvec];
		elseif(cellProp.isInterneuronCell==1)
			fsCount=fsCount+1;
			%isi_FS(fsCount)=nanmean(1./diff(cellProp.spikeTimes));
			isi_FS=[isi_FS ISIvec];
		else
			grayCount=grayCount+1;
			%isi_Gray(grayCount)=nanmean(1./diff(cellProp.spikeTimes));
			%isi_Gray(grayCount)=cellProp.meanRate;
			isi_Gray=[isi_Gray ISIvec];
		end
	end	
	toc
	
	%save cross corrs as either FS, RS, or gray
	figure
	%binWidth=0.1;
	binWidth=0.05;
	expBinWidth=0.5;
	expBinWidth=0.025;
	edgeExps=-3:expBinWidth:2;

	edges=10.^edgeExps
	histogram(isi_RS,edges)
	%histogramLog(isi_RS,edges)
	
	set(gca,'xscale','log')
	xlabel('ISI (sec)')
	ylabel('ISI count')

	%title('RS avg inst. frequency')
	title([groupName '-RS ISI distribution'])
	%title(sprintf('%s-RS avg firing rate: %.2f +/- %.2f',groupName,nanmean(isi_RS),nanstd(isi_RS)))
	
	set(gca,'FontSize',16)	
	figure
	histogram(isi_FS,edges)
	%histogramLog(isi_FS,edges)
	set(gca,'xscale','log')
	xlabel('ISI (sec)')
	ylabel('ISI count')
	title([groupName '-FS ISI distribution'])
	
	%title('FS avg inst. frequency')
	%title(sprintf('%s-FS avg firing rate: %.2f +/- %.2f',groupName,nanmean(isi_FS),nanstd(isi_FS)))
	%save avg cross cors for FS, RS, gray
	set(gca,'FontSize',16)	
	disp(sprintf('FS avg.: %.2f Hz',nanmean(isi_FS)))
	disp(sprintf('RS avg.: %.2f Hz',nanmean(isi_RS)))

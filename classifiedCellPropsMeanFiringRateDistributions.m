function [] = classifiedCellPropsDynamics(cellPropDir,groupName)
	%make sure input directory mat files are all same group
	%cellPropFilePaths=getFilePathsRegex(cellPropDir,'*_seizure1_*.mat')
	cellPropFilePaths=getFilePathsRegex(cellPropDir,'*.mat')
	%cellPropFilePaths=getFilePathsRegex(cellPropDir,'*MG49*_seizure1_*.mat')
	%cellPropFilePaths=getFilePathsRegex(cellPropDir,'*.mat')

	rsCount=0;
	fsCount=0;
	grayCount=0;
	clear avgInstFR_RS	
	clear avgInstFR_FS
	clear avgInstFR_Gray
	
	tic
	disp('getting firing rates by cell type.....')
	%save firing rate distributions by cell type	
	for cellIdx=1:length(cellPropFilePaths)
		cellProp=load(cellPropFilePaths{cellIdx});
		if(cellProp.isInterneuronCell==0)
			rsCount=rsCount+1;
			%avgInstFR_RS(rsCount)=nanmean(1./diff(cellProp.spikeTimes));
			avgInstFR_RS(rsCount)=cellProp.meanRate;
		elseif(cellProp.isInterneuronCell==1)
			fsCount=fsCount+1;
			%avgInstFR_FS(fsCount)=nanmean(1./diff(cellProp.spikeTimes));
			avgInstFR_FS(fsCount)=cellProp.meanRate;
		else
			grayCount=grayCount+1;
			%avgInstFR_Gray(grayCount)=nanmean(1./diff(cellProp.spikeTimes));
			avgInstFR_Gray(grayCount)=cellProp.meanRate;
		end
	end	
	toc
	
	%save cross corrs as either FS, RS, or gray
	figure
	%binWidth=0.1;
	binWidth=0.05;
	edges=(0:binWidth:10);
	expBinWidth=0.5;
	expBinWidth=0.3;
	edgeExps=-3:expBinWidth:2;

	edges=10.^edgeExps
	histogram(avgInstFR_RS,edges)
	%histogramLog(avgInstFR_RS,edges)
	
	set(gca,'xscale','log')
	xlabel('Avg Firing Rate (Hz)')
	ylabel('Cell count')
	%title('RS avg inst. frequency')
	title('RS avg firing rate')
	title(sprintf('%s-RS avg firing rate: %.2f +/- %.2f',groupName,nanmean(avgInstFR_RS),nanstd(avgInstFR_RS)))
	
	set(gca,'FontSize',16)	
	figure
	histogram(avgInstFR_FS,edges)
	%histogramLog(avgInstFR_FS,edges)
	set(gca,'xscale','log')
	xlabel('Avg Firing Rate (Hz)')
	ylabel('Cell count')
	%title('FS avg inst. frequency')
	title(sprintf('%s-FS avg firing rate: %.2f +/- %.2f',groupName,nanmean(avgInstFR_FS),nanstd(avgInstFR_FS)))
	%save avg cross cors for FS, RS, gray
	set(gca,'FontSize',16)	
	disp(sprintf('FS avg.: %.2f Hz',nanmean(avgInstFR_FS)))
	disp(sprintf('RS avg.: %.2f Hz',nanmean(avgInstFR_RS)))

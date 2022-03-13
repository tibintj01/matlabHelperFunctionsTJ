szIDs={'MG49-seizure45','MG49_seizure36', 'MG49-seizure43', 'MG63_seizure1-4','BW9_seizure1-3','RIHE1'}

for szIdx=1:6
	szID=szIDs{szIdx};
	%cellPropPaths=getRegexFilePaths('/nfs/turbo/lsa-ojahmed/compiledDir_01-Nov-2017_23-47-00_getCellWaveformPropsOmar_Tibin/seizureCellWaveformProperties/',sprintf('*cell_properties*%s*.mat',szID));
	%cellDir='/nfs/turbo/lsa-ojahmed/compiledDir_01-Nov-2017_23-47-00_getCellWaveformPropsOmar_Tibin/seizureCellWaveformProperties/';
	cellDir='/home/tibintj/pbsTest/compiledDir_17-Nov-2017_13-07-44_getCellWaveformPropsOmar_Tibin/seizureCellWaveformProperties-RawWaveforms';
	cellPropPaths=getRegexFilePaths(cellDir,sprintf('*cell_properties*%s*.mat',szID));

	maxDescResiduals=NaN(length(cellPropPaths),1);
	avgSpikeTPwidth=NaN(length(cellPropPaths),1);
	avgOvershoot=NaN(length(cellPropPaths),1);
	avgOvershootTrough=NaN(length(cellPropPaths),1);
	avgOvershootInd=NaN(length(cellPropPaths),1);
	avgOvershootTroughInd=NaN(length(cellPropPaths),1);

	
	[dummy]=getAvgWaveformAcrossCells(cellDir,szID);

	for cellIdx=1:length(cellPropPaths)
		cellIdx/length(cellPropPaths)

		badCellList=getBadCellList(cellPropPaths{cellIdx});
		fileName=getFileNameFromPath(cellPropPaths{cellIdx});
		cellName=getStrUpTo(fileName,'_');

		if(max(ismember(badCellList,cellName))==0)
			[maxDescResiduals(cellIdx),avgSpikeTPwidth(cellIdx),avgOvershoot(cellIdx),avgOvershootTrough(cellIdx),avgOvershootInd(cellIdx),avgOvershootTroughInd(cellIdx)]=plotSingleSpikesForFile(cellPropPaths{cellIdx},szID);	
		end
	end

	close all
	plot(maxDescResiduals,avgSpikeTPwidth,'ko', 'MarkerSize',1)
	title([szID])

	xlabel('Max desc. diff from pop. avg. (frac. of peak)')
	ylabel('Trough-to-peak width (bins)')

	saveas(gcf,sprintf('MaxResidualTPspace-%s.tif',szID))

	close all
	plot(avgOvershoot,avgSpikeTPwidth,'ko', 'MarkerSize',1)
	title([szID])

	xlabel('Mean diff from pop. avg. (frac. of peak)')
	ylabel('Trough-to-peak width (bins)')

	saveas(gcf,sprintf('MeanResidualTPspace-%s.tif',szID))

	close all
	plot(avgOvershoot,-avgOvershootTrough,'ko', 'MarkerSize',1)
	title([szID])

	hold on
	plot([min(avgOvershoot)*1.2 max(avgOvershoot)],[0 0],'r--')
	plot([0 0],[min(-avgOvershootTrough)*1.2 max(-avgOvershootTrough)],'r--')

	xlim([min(avgOvershoot)*1.2 max(avgOvershoot)])
	ylim([min(-avgOvershootTrough)*1.2 max(-avgOvershootTrough)])


	ylabel('-Mean asc. diff from pop. avg. (frac. of peak)')
	xlabel('Mean desc. diff from pop. avg. (frac. of trough)')
	daspect([ 1 1 1])

	saveas(gcf,sprintf('MeanResidualsAscAndDescSpace-%s.tif',szID))
	
	close all
	plot(avgOvershootInd,-avgOvershootTroughInd,'ko', 'MarkerSize',1)
	title([szID])

	hold on
	plot([min(avgOvershoot)*1.2 max(avgOvershoot)],[0 0],'r--')
	plot([0 0],[min(-avgOvershootTrough)*1.2 max(-avgOvershootTrough)],'r--')

	xlim([min(avgOvershoot)*1.2 max(avgOvershoot)])
	ylim([min(-avgOvershootTrough)*1.2 max(-avgOvershootTrough)])


	ylabel('-Mean asc. diff from pop. avg. (frac. of peak)')
	xlabel('Mean desc. diff from pop. avg. (frac. of trough)')
	daspect([ 1 1 1])

	saveas(gcf,sprintf('MeanResidualsAscAndDescIndSpace-%s.tif',szID))

end


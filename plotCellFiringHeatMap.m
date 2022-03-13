function [] = plotCellFiringHeatMap(cells,binWidth,timeLims)

	if(~exist('binWidth'))
		binWidth=1;
	end
	if(~exist('timeLims'))
		minTime=0;
		maxTime=cells{1}.tend;
	else
		minTime=timeLims(1);
		maxTime=timeLims(2);
	end

	firingRateMatrix=zeros(size(length(cells),ceil((maxTime-minTime)/binWidth)));
        for cellIdx=1:length(cells)
                cell=cells{cellIdx};

                if(cell.qualityRating<3)
                        continue
                end

		spikeTimes=cell.spikeTimes;


		[binCenters,firingRate]=getDynamicFiringRate(spikeTimes,binWidth,minTime,maxTime); 	

		firingRateMatrix(cellIdx,1:length(firingRate))=firingRate(:)';
		
	end
	meanFiringRates=mean(firingRateMatrix,2);
	[dummyVar,sortedMeanIdxes]=sort(meanFiringRates,'descend');
	firingRateMatrix=firingRateMatrix(sortedMeanIdxes,:);

	figure
	cmin=prctile(firingRateMatrix(:),5);
	cmax=prctile(firingRateMatrix(:),99.99);
	imagesc(firingRateMatrix,[cmin cmax])
	cbarH=colorbar
	ylabel(cbarH,'Firing rate (Hz)')

	numLabels=10;
	labelIdxes=round(linspace(1,length(binCenters),numLabels));
	secToMin=1/60;
	xticklabels=binCenters(labelIdxes);
	xticks=linspace(1,size(firingRateMatrix,2),length(xticklabels))
	set(gca,'XTick',xticks,'XTickLabel',xticklabels)

	%currAx=gca;
	%currAx.XAxis.TickLabelFormat='%.1f';

	xlabel('Time from start (sec)')
	ylabel('Sorted cell no.')

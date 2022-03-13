function [descendingPhaseMaxResidual,avgWaveformTPwidth,avgOvershoot,avgOvershootTrough,avgOvershootsInd,avgOvershootsIndTrough] = plotSingleSpikesForFile(filePath,szID)
splineDt=0.1;
%cellPropI=load('/nfs/turbo/lsa-ojahmed/compiledDir_01-Nov-2017_23-47-00_getCellWaveformPropsOmar_Tibin/seizureCellWaveformProperties/9c_cell_properties_MG49-seizure45.mat')
%cellPropI=load('/nfs/turbo/lsa-ojahmed/compiledDir_01-Nov-2017_23-47-00_getCellWaveformPropsOmar_Tibin/seizureCellWaveformProperties/79a_cell_properties_MG49-seizure45.mat')
%cellPropE=load('/nfs/turbo/lsa-ojahmed/compiledDir_01-Nov-2017_23-47-00_getCellWaveformPropsOmar_Tibin/seizureCellWaveformProperties/9b_cell_properties_MG49-seizure45.mat')
%cellPropE=load('/nfs/turbo/lsa-ojahmed/compiledDir_01-Nov-2017_23-47-00_getCellWaveformPropsOmar_Tibin/seizureCellWaveformProperties/64b_cell_properties_MG49-seizure45.mat')

%eFileName='62a_cell_properties_MG49-seizure45.mat';

%iFileName='22c_cell_properties_MG49-seizure45.mat';
%eFileName='22b_cell_properties_MG49-seizure45.mat';
%iFileName='29a_cell_properties_MG49-seizure45.mat';
%eFileName='29b_cell_properties_MG49-seizure45.mat';
%iFileName='68a_cell_properties_MG49-seizure45.mat';
%eFileName='62a_cell_properties_MG49-seizure45.mat';
%iFileName='79a_cell_properties_MG49-seizure45.mat';
%eFileName='82a_cell_properties_MG49-seizure45.mat';
%iFileName='9c_cell_properties_MG49-seizure45.mat';
%eFileName='6a_cell_properties_MG49-seizure45.mat';
%iFileName='34a_cell_properties_MG49-seizure45.mat';
%eFileName='36a_cell_properties_MG49-seizure45.mat';

%eCellID=['e' eFileName(1:3)];
%iCellID=['i' iFileName(1:3)];
fileName=getFileNameFromPath(filePath);
cellID=fileName(1:3);

%cellPropI=load(sprintf('/nfs/turbo/lsa-ojahmed/compiledDir_01-Nov-2017_23-47-00_getCellWaveformPropsOmar_Tibin/seizureCellWaveformProperties/%s',iFileName));
%cellPropE=load(sprintf('/nfs/turbo/lsa-ojahmed/compiledDir_01-Nov-2017_23-47-00_getCellWaveformPropsOmar_Tibin/seizureCellWaveformProperties/%s',eFileName));

cellProp=load(filePath);


waveforms=getGoodWaveforms(filePath);

numSpikesMin=3;
if(size(waveforms,2)<numSpikesMin)
	descendingPhaseMaxResidual=NaN;
	avgWaveformTPwidth=NaN;
	avgOvershoot=NaN;
	avgOvershootTrough=NaN;
	avgOvershootsInd=NaN;
	avgOvershootsIndTrough=NaN;
	disp(sprintf('Less than %d good spikes! Skipping....',numSpikesMin))
	return
end

%allWaveforms=[iWaveforms eWaveforms];

%allAvgWaveform=mean(allWaveforms,2);
%cellAvgSpikeFileName='allCellAvgWaveformMG49-seizure45.mat';
cellAvgSpikeFileName=sprintf('allCellAvgWaveform%s.mat',szID);
%cellDir='/nfs/turbo/lsa-ojahmed/compiledDir_01-Nov-2017_23-47-00_getCellWaveformPropsOmar_Tibin/seizureCellWaveformProperties/';
cellDir=getDirNameFromPath(filePath);
[dummy]=getAvgWaveformAcrossCells(cellDir,szID);

%allAvgWaveformFile=load(sprintf('/nfs/turbo/lsa-ojahmed/compiledDir_01-Nov-2017_23-47-00_getCellWaveformPropsOmar_Tibin/seizureCellWaveformProperties/%s',cellAvgSpikeFileName));
%allAvgWaveformFile=load(sprintf(cellDir,cellAvgSpikeFileName));
allAvgWaveformFile=load(fullfile(cellDir,cellAvgSpikeFileName));

allAvgWaveform=allAvgWaveformFile.allCellAvgWaveform;

allAvgWaveformTrough=minNormalize(allAvgWaveform);

allAvgWaveform=rangeNormalize(allAvgWaveform);
waveformsTrough=minNormalizeCols(waveforms);

waveforms=rangeNormalizeCols(waveforms);
%eWaveforms=rangeNormalizeCols(eWaveforms);

%waveformsRes=waveforms-repmat(allAvgWaveform,1,size(waveforms,2));

waveformsTroughAligned=zeros(size(waveforms));
waveformsTroughAlignedPnorm=zeros(size(waveforms));
for spikeIdx=1:size(waveforms,2)
	waveformsTroughAligned(:,spikeIdx)=troughAlignTo2nd(waveformsTrough(:,spikeIdx),allAvgWaveformTrough);
	waveformsTroughAlignedPnorm(:,spikeIdx)=troughAlignTo2nd(waveforms(:,spikeIdx),allAvgWaveform);
end
waveformsResTrough=waveformsTroughAligned-repmat(allAvgWaveformTrough,1,size(waveforms,2));

waveformsRes=waveformsTroughAlignedPnorm-repmat(allAvgWaveform,1,size(waveforms,2));

close all
startIdx=1;
%numSpikes=60;

%average residue for each spike...
%wResIndividual=mean(waveformsRes,2);
%wResIndividualTrough=mean(waveformsResTrough,2);

%or get residue of average spike?
normAvgWaveform=rangeNormalize(mean(waveforms,2));
normAvgWaveformTrough=minNormalize(mean(waveforms,2));

[avgMinCell,avgMinCellIdx]=min(normAvgWaveform);

%align based on descending phase of spike with average waveform
normAvgWaveform=xcAlign(normAvgWaveform,allAvgWaveform,avgMinCellIdx);
%fds



[avgMinCell,avgMinCellBin]=min(normAvgWaveform);
[avgMinCellTrough,avgMinCellBinTrough]=min(normAvgWaveformTrough);
[maxAvgVal,maxAvgBin]=max(normAvgWaveform);

[avgMinPop,avgMinPopBin]=min(allAvgWaveform);
[avgMinPopTrough,avgMinPopBinTrough]=min(allAvgWaveformTrough);
[maxAvgPop,avgMaxPopBin]=max(allAvgWaveform);

plotTPcomparison=0;
%if(length(findstr(cellID,'14a'))>0)
%	plotTPcomparison=1;
%end
if(plotTPcomparison)
	avgTP=maxAvgBin-avgMinCellBin;
	avgTPoriginalBin=splineToOriginalBin(avgTP);

	figure; histogram(cellProp.spikeTPWidth)
	title(sprintf('Avg TP bin: %d',avgTPoriginalBin))
	saveas(gcf,sprintf('avgWaveformTPvsAvgTP%s-%s.tif',cellID,szID))
%	fds
end

%try aligned peak to see descending shape difference independent of timing
%peakAlignShift=(avgMaxPopBin-maxAvgBin);
%normAvgWaveform=shiftVals(normAvgWaveform,peakAlignShift);
[normAvgWaveform,peakAlignShift]=xcAlignUpTo(normAvgWaveform,allAvgWaveform,maxAvgBin);


waveformsResXCaligned=zeros(size(waveformsRes));
for spikeIdx=1:size(waveformsRes,2)
	[waveformsResXCaligned(:,spikeIdx),peakAlignShift]=xcAlignUpTo(waveformsRes(:,spikeIdx),allAvgWaveform,avgMaxPopBin);
	[dummy,indSpikeMaxIdx(spikeIdx)]=max(waveformsResXCaligned(:,spikeIdx));
end


%[normAvgWaveformTrough,peakAlignShiftTrough]=xcAlignUpTo(normAvgWaveformTrough,allAvgWaveformTrough,maxAvgBin);
peakAlignShiftTrough=(avgMinPopBinTrough-avgMinCellBinTrough);
normAvgWaveformTrough=shiftVals(normAvgWaveformTrough,peakAlignShiftTrough);

wRes=normAvgWaveform(:)-allAvgWaveform(:);
wResTrough=normAvgWaveformTrough(:)-allAvgWaveformTrough(:);

%[leftIdxI,leftExtremeI,rightIdxI,rightExtremeI]=getLeftRightExtrema(iRes);
%[leftIdxE,leftExtremeE,rightIdxE,rightExtremeE]=getLeftRightExtrema(eRes);
%maxOriginalBin=nanmean(cellProp.spikeMaxIndex);
%maxHighBin=round(maxOriginalBin/splineDt-1+1/splineDt); %convert back to splined resolution

%binThresh=min(length(wRes),maxAvgBin+10);
binThresh=min(length(wRes),avgMaxPopBin+10);
binThreshTrough=min(length(wResTrough),avgMinPopBin+10);

popWaveBetween=allAvgWaveformTrough(avgMinPopBin:avgMaxPopBin);

[dummyZero,popZeroCrossBinFromMin]=min(abs(popWaveBetween));
popZeroCrossBin=popZeroCrossBinFromMin+avgMinPopBin-1;

[rightMinVal,rightMinIdx]=getRightMin(wRes,binThresh);
if(peakAlignShift>=0)
	avgOvershoot=nanmean(wRes(binThresh:end));
else
	avgOvershoot=nanmean(wRes(binThresh:(end+peakAlignShift)));
end
	avgOvershootTrough=nanmean(wResTrough(binThreshTrough:popZeroCrossBin));

%avgOvershootsInd=nanmean(nanmean(waveformsRes(binThresh:end,:),1));
abovePerc=0;
bottomPerc=30;
avgOvershootsInd=abovePercMean(nanmean(waveformsRes(binThresh:end,:),1),abovePerc);
avgOvershootsIndTrough=-abovePercMean(-nanmean(waveformsResTrough(binThreshTrough:popZeroCrossBin,:),1),abovePerc);
%avgOvershootsInd=bottomPercMean(nanmean(waveformsRes(binThresh:end,:),1),bottomPerc);

%avgOvershootsIndTrough=nanmean(nanmean(waveformsResTrough(binThreshTrough:popZeroCrossBin,:),1));
%avgOvershootsInd=nanmean(nanmean(waveformsResXCaligned(binThresh:end,:),1));

%avgOvershootsIndAll=zeros(length(indSpikeMaxIdx),1);
%for i=1:length(indSpikeMaxIdx)
%	avgOvershootsIndAll(i)=(nanmean(waveformsResXCaligned((indSpikeMaxIdx(i)+10):end,i),1));
%end

%avgOvershootsInd=nanmean(avgOvershootsIndAll);

%avgOvershootsIndTrough=nanmean(nanmean(waveformsResTrough(binThreshTrough:popZeroCrossBin,:),1));
%avgOvershootsIndAll=zeros(length(indSpikeMaxIdx),1);

%avgOvershootsIndTrough=nanmean(nanmean(waveformsResTrough(binThreshTrough:popZeroCrossBin,:),1));
%avgOvershootsIndTrough=-bottomPercMean(-nanmean(waveformsResTrough(binThreshTrough:popZeroCrossBin,:),1),bottomPerc);

%avgOvershootsIndTroughs=zeros(size(indSpikeMaxIdx));

%for i=1:length(indSpikeMaxIdx)
	%avgOvershootsIndAll(i)=(nanmean(waveformsResXCaligned((indSpikeMaxIdx(i)+10):end,i),1));
%	avgOvershootsIndTroughs(i)=(nanmean(waveformsResTrough(binThreshTrough:indSpikeMaxIdx(i),i),1));
%end
%avgOvershootsIndTrough=nanmean(avgOvershootsIndTroughs);

descendingPhaseMaxResidual=rightMinVal;
avgWaveformTPwidth=(splineToOriginalBin(maxAvgBin)-splineToOriginalBin(avgMinCellBin));
if(avgWaveformTPwidth<0)
	descendingPhaseMaxResidual=NaN;
	avgWaveformTPwidth=NaN;
	avgOvershoot=NaN;
	avgOvershootTrough=NaN;
	avgOvershootsInd=NaN;
	avgOvershootsIndTrough=NaN;
	disp('TP less than zero! Skipping....')
	return
end
save(sprintf('avgResidualFeatures%s_%s.mat',cellID,szID),'descendingPhaseMaxResidual','avgWaveformTPwidth','avgOvershoot','avgOvershootTrough','avgOvershootsInd','avgOvershootsIndTrough')

subplot(1,2,1)
%plot(allAvgWaveform)
maxSpikeIdx=min(size(waveforms,2),(startIdx+4));
plot(normAvgWaveformTrough,'b')
hold on
plot(wResTrough,'r')
plot(allAvgWaveformTrough,'k')

%plot(waveforms(:,startIdx:maxSpikeIdx))
plot(popZeroCrossBin,avgOvershootTrough,'go')
%plot(popZeroCrossBin,avgOvershootsIndTrough,'bo')
plot(round(mean(indSpikeMaxIdx)),avgOvershootsIndTrough,'bo')
%title(sprintf('Cell %s spike examples',cellID))
title(sprintf('Cell %s spike diff; trough normalized',cellID))

addCustomLegendTibin({'go'},{'Avg. Diff. Trough to 0-crossing'})
yMax=1.1;
ylim([-yMax yMax])
subplot(1,2,2)
%plot(mean(iWaveforms(:,startIdx:(startIdx+numSpikes-1)),2))
%plot(mean(waveforms,2),'b')
plot(normAvgWaveform,'b')
%plot(iWaveformsRes(:,startIdx:(startIdx+numSpikes-1)))

hold on
plot(wRes,'r')
plot(allAvgWaveform,'k')
%plot(leftIdxI,leftExtremeI,'bo')
%plot(rightIdxI,rightExtremeI,'ro')
%plot(rightMinIdx,rightMinVal,'ro')
plot(binThresh,avgOvershoot,'go')
plot(binThresh,avgOvershootsInd,'bo')

title(sprintf('Cell %s spike diff; peak normalized',cellID))
%plot(iWaveformsRes(:,startIdx:(startIdx+numSpikes-1)))
%title('I cell spike residues')

%ylim([-40 50])

ylim([-yMax yMax])
%addCustomLegendTibin({'b-','k-','r-','ro'},{'Avg. Waveform Cell','Avg. Waveform Pop.', 'Residual','2nd half min'})
addCustomLegendTibin({'k-','b-','r-','go'},{'Pop. Avg. Waveform','Avg. Waveform Cell', 'Difference','Avg. Diff. Peak to end'})

%plot(goodWaveformsInterp(:,1:numSpikes),'ko','MarkerSize',1)

%saveas(gcf,sprintf('singleSpikePlots%sAND%s.tif',iCellID,eCellID))
%saveas(gcf,sprintf('singleSpikePlotsResiduals%s_%s.tif',cellID,szID))
saveas(gcf,sprintf('peakAndTroughAlignedResiduals%s_%s.tif',cellID,szID))

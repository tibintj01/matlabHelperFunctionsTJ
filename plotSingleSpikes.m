%cellPropI=load('/nfs/turbo/lsa-ojahmed/compiledDir_01-Nov-2017_23-47-00_getCellWaveformPropsOmar_Tibin/seizureCellWaveformProperties/9c_cell_properties_MG49-seizure45.mat')
%cellPropI=load('/nfs/turbo/lsa-ojahmed/compiledDir_01-Nov-2017_23-47-00_getCellWaveformPropsOmar_Tibin/seizureCellWaveformProperties/79a_cell_properties_MG49-seizure45.mat')
%cellPropE=load('/nfs/turbo/lsa-ojahmed/compiledDir_01-Nov-2017_23-47-00_getCellWaveformPropsOmar_Tibin/seizureCellWaveformProperties/9b_cell_properties_MG49-seizure45.mat')
%cellPropE=load('/nfs/turbo/lsa-ojahmed/compiledDir_01-Nov-2017_23-47-00_getCellWaveformPropsOmar_Tibin/seizureCellWaveformProperties/64b_cell_properties_MG49-seizure45.mat')

%eFileName='62a_cell_properties_MG49-seizure45.mat';

%iFileName='22c_cell_properties_MG49-seizure45.mat';
%eFileName='22b_cell_properties_MG49-seizure45.mat';
iFileName='29a_cell_properties_MG49-seizure45.mat';
eFileName='29b_cell_properties_MG49-seizure45.mat';
iFileName='68a_cell_properties_MG49-seizure45.mat';
eFileName='62a_cell_properties_MG49-seizure45.mat';
iFileName='79a_cell_properties_MG49-seizure45.mat';
eFileName='82a_cell_properties_MG49-seizure45.mat';
iFileName='9c_cell_properties_MG49-seizure45.mat';
eFileName='6a_cell_properties_MG49-seizure45.mat';
iFileName='34a_cell_properties_MG49-seizure45.mat';
eFileName='36a_cell_properties_MG49-seizure45.mat';

eCellID=['e' eFileName(1:3)];
iCellID=['i' iFileName(1:3)];

cellPropI=load(sprintf('/nfs/turbo/lsa-ojahmed/compiledDir_01-Nov-2017_23-47-00_getCellWaveformPropsOmar_Tibin/seizureCellWaveformProperties/%s',iFileName));
cellPropE=load(sprintf('/nfs/turbo/lsa-ojahmed/compiledDir_01-Nov-2017_23-47-00_getCellWaveformPropsOmar_Tibin/seizureCellWaveformProperties/%s',eFileName));

iWaveforms=getGoodWaveforms(cellPropI,iFileName);
eWaveforms=getGoodWaveforms(cellPropE,eFileName);

%allWaveforms=[iWaveforms eWaveforms];

%allAvgWaveform=mean(allWaveforms,2);
cellAvgSpikeFileName='allCellAvgWaveformMG49-seizure45.mat';

allAvgWaveformFile=load(sprintf('/nfs/turbo/lsa-ojahmed/compiledDir_01-Nov-2017_23-47-00_getCellWaveformPropsOmar_Tibin/seizureCellWaveformProperties/%s',cellAvgSpikeFileName));

allAvgWaveform=allAvgWaveformFile.allCellAvgWaveform;

allAvgWaveform=rangeNormalize(allAvgWaveform);
iWaveforms=rangeNormalizeCols(iWaveforms);
eWaveforms=rangeNormalizeCols(eWaveforms);

iWaveformsRes=iWaveforms-repmat(allAvgWaveform,1,size(iWaveforms,2));
eWaveformsRes=eWaveforms-repmat(allAvgWaveform,1,size(eWaveforms,2));
close all
startIdx=10;
numSpikes=60;

iRes=mean(iWaveformsRes,2);
eRes=mean(eWaveformsRes,2);

%[leftIdxI,leftExtremeI,rightIdxI,rightExtremeI]=getLeftRightExtrema(iRes);
%[leftIdxE,leftExtremeE,rightIdxE,rightExtremeE]=getLeftRightExtrema(eRes);

binThresh=200;
[rightMinVal_I,rightMinIdx_I]=getRightMin(iRes,binThresh);
[rightMinVal_E,rightMinIdx_E]=getRightMin(eRes,binThresh);

subplot(1,3,1)
plot(allAvgWaveform)
title('All avg waveform')

subplot(1,3,2)
%plot(mean(iWaveforms(:,startIdx:(startIdx+numSpikes-1)),2))
plot(mean(iWaveforms,2),'b')
%plot(iWaveformsRes(:,startIdx:(startIdx+numSpikes-1)))

hold on
plot(iRes,'r')
plot(allAvgWaveform,'k')
%plot(leftIdxI,leftExtremeI,'bo')
%plot(rightIdxI,rightExtremeI,'ro')
plot(rightMinIdx_I,rightMinVal_I,'ro')

title('I cell spike avg residuals')
%plot(iWaveformsRes(:,startIdx:(startIdx+numSpikes-1)))
%title('I cell spike residues')

%ylim([-40 50])
ylim([-0.6 0.6])

subplot(1,3,3)
%plot(mean(eWaveforms(:,startIdx:(startIdx+numSpikes-1)),2))
plot(mean(eWaveforms,2),'b')
hold on
%plot(eWaveformsRes(:,startIdx:(startIdx+numSpikes-1)))
%title('E cell spike residues')

plot(eRes,'r')
plot(allAvgWaveform,'k')
%plot(leftIdxE,leftExtremeE,'bo')
%plot(rightIdxE,rightExtremeE,'ro')
plot(rightMinIdx_E,rightMinVal_E,'ro')
title('E cell spike avg residues')
%ylim([-40 50])
ylim([-0.6 0.6])

%addCustomLegendTibin({'b-','k-','r-','ro'},{'Avg. Waveform Cell','Avg. Waveform Pop.', 'Residual','2nd half min'})
addCustomLegendTibin({'b-','k-','r-'},{'Avg. Waveform Cell','Avg. Waveform Pop.', 'Residual'})

%plot(goodWaveformsInterp(:,1:numSpikes),'ko','MarkerSize',1)

saveas(gcf,sprintf('singleSpikePlots%sAND%s.tif',iCellID,eCellID))

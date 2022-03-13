%Exploratory script to investigate spiking and LFP in pt MG49
close all

%cells=getCellWaveformTimingMatrixForPt();
%cellSaveName='sess3CellsWithExtendedRawWaveforms.mat';
%newCells=addRawWaveforms(cells,cellSaveName);
%[rawWavforms,newCell]=getExtendedRawSpikeWaveforms(cells{2});
%plot(newCell.rawSpikeWaveforms')
%fds

cellSaveName='sess3CellsWithExtendedRawWaveforms.mat';
cellData=load(cellSaveName);
cells=cellData.newCells;

%cellData=matfile('Sess3CellsWithRawWaveforms.mat');
%fds
%cells=cellData.newCells;

%correlate filtered LFP and cell waveform amplitude?
exCell=cells{6}
disp('filtering lfp.......')
tic
%[lfpDecFact,lfpForCell]=getLFPforCell(exCell);
[lfpForCell,lfpFs]=getDecLFPforCell(exCell);
disp('done')
toc
lfpDecFact=round(exCell.Fs/lfpFs);
spikeIdxes=exCell.spikeTimes*exCell.Fs;
lfpAtSpikes=lfpForCell(round(spikeIdxes/lfpDecFact));

%spikeAmps=max(exCell.waveForms,[],2);
spikeAmps=prctile(exCell.waveForms,90,2);
ISIs=diff(exCell.spikeTimes);

figure
disp('plotting.....')
subplot(2,1,1)
plot(lfpAtSpikes,spikeAmps,'ko','MarkerSize',1)
xlabel('Same-Channel LFP Amp')
ylabel('10th percentile Spike Amp')

subplot(2,1,2)
plot(lfpAtSpikes(2:end),ISIs,'ko','MarkerSize',1)
xlabel('Same-Channel LFP Amp')
ylabel('ISI')


disp('done')
%newCells=addRawWaveforms(cells);

%exCell=cells{1}
%find raw waveforms from raw lfp based on spike times
%[rawWaveforms,newCell]=getRawSpikeWaveforms(exCell);
%plot(rawWaveforms)

%compute clustering on raw and filtered

%eiClustering(cells,0)
%eiClustering(cells,1)
%eiClustering(newCells)

%timeLims=[60 120];
%plotCellFiringHeatMap(cells,0.005,timeLims)

%timeLims=[60 120];
%plotCellFiringHeatMap({exCell},0.005,timeLims)

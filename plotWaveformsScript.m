cellSaveName='sess3CellsWithExtendedRawWaveforms.mat';
cellData=load(cellSaveName);
cells=cellData.newCells;

useRaw=0;
%useRaw=1;
%useRaw=2;
plotWaveformsByCell(cells,useRaw)


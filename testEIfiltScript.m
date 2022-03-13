cellSaveName='sess3CellsWithExtendedRawWaveforms.mat';
cellData=load(cellSaveName);
cells=cellData.newCells;

eiClustering(cells,1)

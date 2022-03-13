
%read in cells
cellSaveName='sess3CellsWithExtendedRawWaveforms.mat';
cellData=load(cellSaveName);
cells=cellData.newCells;


%inhibitory cells: 108, 10

cellIdxes=[108 10 36 66];
phenotypes={'FS', 'FS', 'RS','RS'}
maxNumSpike=100000;

for i=1:length(cellIdxes)
%for i=2:2
	cellIdx=cellIdxes(i);
	phenotype=phenotypes{i};

	cell=cells{cellIdx}
	titleStr=sprintf('Cell-%d-%s',cellIdx,phenotype)
	plotPeriSpikeForCell(cell,titleStr,maxNumSpike)
end

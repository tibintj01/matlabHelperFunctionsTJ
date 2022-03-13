ptDir='MG49';
sessionNum=3;
ch=6;

cycleSegmentedSpikeProperties=load('singleCyclebroadbandAlphaDurationNormed-Ch06SpikeProps.mat');
cycleSegmentedSpikeProperties=cycleSegmentedSpikeProperties.cycleSegmentedSpikeProperties;
%getCycleSegmentedSpikeProperties(ptDir,sessionNum,ch,'broadbandAlpha');

cycleWaveformSubspace=load('singleCycleBroadbandAlpha-DurationNormed-AllTroughs-Ch06SVD.mat');

durationData=load('cycleMatrixBroadbandZscore-MG49-Alpha-zScoredCycles-DurationNormalized-Ch06.mat');

cycleDurations=durationData.actualDurations;

%each cycle has a sigma(i)*v(i) amount of u(i)
U3=cycleWaveformSubspace.U(:,1:3);
V3=cycleWaveformSubspace.V(:,1:3);
sigmas=diag(cycleWaveformSubspace.S(:,1:3));

numCycles=size(V3,1)

figure

cycleSpikeCounts=cycleSegmentedSpikeProperties.cycleSpikeCounts;
numCells=size(cycleSpikeCounts,1);


cycleHasSpike=cycleSpikeCounts>0;

cycleComponentStrengthBin=zeros(3,numCycles);

numBinsPerComp=10;
for i=1:3
	compStrengthVector=sigmas(i)*V3(:,i);
	binEdges=linspace(prctile(compStrengthVector,5),prctile(compStrengthVector,95),numBinsPerComp);
	cycleComponentStrengthBin(i,:)=discretize(compStrengthVector,binEdges);
end

%spikeProbInBin=zeros(numBinsPerComp,numBinsPerComp,numBinsPerComp)
spikeProbInBin=zeros(numBinsPerComp,numBinsPerComp);
%spikeProbInBin=zeros(numBinsPerComp,1);
%for j=1:3
meanCycleDurationInBin=zeros(numBinsPerComp,numBinsPerComp);

for cellNum=1:numCells
	for i=1:numBinsPerComp
		for j=1:numBinsPerComp
			%numCyclesWithSpikeInBin=sum(cycleHasSpike(j,:) & cycleComponentStrengthBin(j,:) == i);
			numCyclesWithSpikeInBin=sum(cycleHasSpike(cellNum,:) & cycleComponentStrengthBin(1,:) == i & cycleComponentStrengthBin(2,:) == j);
			%totalNumCyclesInBin=sum(cycleComponentStrengthBin(j,:) == i);
			totalNumCyclesInBin=sum(cycleComponentStrengthBin(1,:) == i & cycleComponentStrengthBin(2,:) == j);
			%spikeProbInBin(i)=numCyclesWithSpikeInBin/totalNumCyclesInBin;
			%spikeProbInBin(i,j)=numCyclesWithSpikeInBin/totalNumCyclesInBin;
			spikeProbInBin(i,j)=spikeProbInBin(i,j)+numCyclesWithSpikeInBin/totalNumCyclesInBin;
			cycleDurationsInBin=cycleDurations(cycleComponentStrengthBin(1,:) == i & cycleComponentStrengthBin(2,:) == j);
			meanCycleDurationInBin(i,j)=mean(cycleDurationsInBin);
		end
	end

end
spikeProbInBin(i,j)=spikeProbInBin(i,j)/numCells;
figure
%plot(spikeProbInBin)
subplot(2,2,2)
imagesc(spikeProbInBin(1:(end-1),1:(end-1)))
xlabel('2nd waveform component strength bin')
ylabel('1st waveform component strength bin')
title('MG49 Ch06 Alpha Waveform SVD subspace and Spiking')
colormap(parula)
c=colorbar
c.Label.String='Probability of at least 1 spike';
%saveas(gcf,sprintf('testSpikeProbForCompBin_cell%d.tif',cellNum))

daspect([1 1 1])

comp1=U3(:,1)*sigmas(1)*mean(V3(:,1));
comp2=U3(:,2)*sigmas(2)*mean(V3(:,2));

subplot(2,2,1)
plot(linspace(0,1,length(U3(:,1))),comp1)
xlabel('Normalized time in cycle')
ylabel('Z-score')
title('Waveform component 1 (avg magnitude across cycles)')

subplot(2,2,4)
plot(linspace(0,1,length(U3(:,1))),comp2)
xlabel('Normalized time in cycle')
ylabel('Z-score')
title('Waveform component 2 (avg magnitude across cycles)')
set(gcf,'units','normalized','outerposition',[0 0 1 1])


subplot(2,2,3)
imagesc(meanCycleDurationInBin(1:(end-1),1:(end-1)))
xlabel('2nd waveform component strength bin')
ylabel('1st waveform component strength bin')
title('MG49 Ch06 Alpha Waveform Subspace Cycle Durations')
colormap(parula)
c=colorbar
c.Label.String='Mean cycle duration (sec)';

daspect([1 1 1])

saveas(gcf,sprintf('testSpikeProbForCompBin_cellAveraged.tif'))

fds
for componentNum=1:3
	
	figure

	for cellNum=1:numCells
		subplot(1,numCells,cellNum)
		plot(sigmas(componentNum)*V3(:,componentNum),cycleSpikeCounts(cellNum,:),'ko','MarkerSize',1)
		ylabel('Spike count in cycle')
		xlabel(sprintf('strength of component %d',componentNum))
		title(sprintf('Cell %d',cellNum))
	end
end



data=load('cycleMatrixBroadbandZscore-MG49-Alpha-zScoredCycles-DurationNormalized.mat')

cycleMatrix=data.cycleMatrix;

numCycles=size(cycleMatrix,2);

numTroughsFilter=zeros(1,size(cycleMatrix,2));

disp('counting troughs....')
for cycleNum=1:numCycles
	currCycle=cycleMatrix(:,cycleNum);

	[pks]=findpeaks(-currCycle);
	numTroughsFilter(cycleNum)=length(pks.loc);
end

singleTroughCycleNums=find(numTroughsFilter==1);

singleTroughCycles=cycleMatrix(:,singleTroughCycleNums);

size(singleTroughCycles)

figure
plot(singleTroughCycles(:,1:100))
cycleMatrix=singleTroughCycles;

save('cycleMatrixBroadbandZscore-MG49-Alpha-zScoredCycles-DurationNormalized-SingleTrough.mat','cycleMatrix')

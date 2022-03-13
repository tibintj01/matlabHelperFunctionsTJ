%freqBandName='Alpha';
freqBandName='Alpha-DurationNormed-AllTroughs-Ch06';
decFs=30000/15;
%if(~exist(sprintf('singleCycle%sSVD.mat',freqBandName)))
if(~exist(sprintf('singleCycleBroadband%sSVD.mat',freqBandName)))
    disp('loading single cycle matrix....')
    tic
    %data=load('/nfs/turbo/lsa-ojahmed/mg49_singleCycleStats_Par/cycleMatrixMaxNormed-MG49-Alpha-zScoredCycles-100.0-ms-Window.mat')
   % data=load('cycleMatrixBroadbandZscore-MG49-Alpha-zScoredCycles-200.0-ms-Window.mat');
    %data=load('cycleMatrixBroadbandZscore-MG49-Alpha-zScoredCycles-150.0-ms-Window.mat');
    %data=load('cycleMatrixBroadbandZscore-MG49-Alpha-zScoredCycles-100.0-ms-Window.mat');
     data=load('/nfs/turbo/lsa-ojahmed/cycleMatrixBroadbandZscore-MG49-Alpha-zScoredCycles-DurationNormalized.mat')
     %data=load('/nfs/turbo/lsa-ojahmed/cycleMatrixBroadbandZscore-MG49-Alpha-zScoredCycles-DurationNormalized-SingleTrough.mat')
  data=load('nfs/turbo/lsa-ojahmed/spikeDynamicsAnalysisTibin/cycleMatrixBroadbandZscore-MG49-Alpha-zScoredCycles-DurationNormalized-Ch06.mat');
    toc

    singleCycles=data.cycleMatrix;

    tic
    disp('computing SVD')
    %maxNumCols=min(size(singleCycles,2),1e6);
    maxNumCols=min(size(singleCycles,2),2e6);
    
    %midTime=0.1;
   % midTime=0.075;
    %startSVDtime=midTime-0.03;
    %endSVDtime=midTime+0.03;

    %svdStartIdx=round(startSVDtime*decFs);
    %svdEndIdx=round(endSVDtime*decFs);

    %[U,S,V]=svd(singleCycles(svdStartIdx:svdEndIdx,1:maxNumCols),'econ');

    %[U,S,V]=svd(singleCycles(:,1:maxNumCols),'econ');
    %[U,S,V]=svd(singleCycles(:,1:maxNumCols));
    [U,S,V]=svd(singleCycles(:,1:end),'econ');
    toc
    numVectorsSave=20;
    U=U(:,1:numVectorsSave);
    S=S(:,1:numVectorsSave);
    V=V(:,1:numVectorsSave);
    save(sprintf('singleCycleBroadband%sSVD.mat',freqBandName),'U','S','V');

else
    %svdData=load(sprintf('singleCycle%sSVD.mat',freqBandName));
    svdData=load(sprintf('singleCycleBroadband%sSVD.mat',freqBandName));
    U=svdData.U;
    S=svdData.S;
    V=svdData.V;
end

sigma=diag(S);

figure
numComponents=5;
subplot(numComponents+1,1,1)


%plot(sigma(1:20),'r*')
plot((sigma(1:20)),'r*')
xlabel('Component number')
ylabel('Singular value')
title('Alpha: Strength of components across cycles')
timeAxis=(1:size(U,1))/decFs;
for i=1:numComponents
   subplot(numComponents+1,2,2*i+1)
   
   plot(timeAxis,U(:,i))
   %xlim([1 length(U(:,i))])
   xlabel('Time in cycle')
   ylabel('zscore LFP')
   title(sprintf('Waveform component %d',i))
   
   subplot(numComponents+1,2,2*i+2)
    midTime=0.1;
   binWidth=0.1;

	%histogram(V(1:1.9e5,i)'*sigma(i))
	histogram(V(1:maxNumCols,i)'*sigma(i))
   %[N,edges]=histcounts(V(:,i)*sigma(i),-0.002:binWidth:0.002)
    %binCenters=edges(2:end)-binWidth/2;
 
	%plot(binCenters,N) 
   xlabel('Component coeff to reproduce data')
   ylabel('Count')

end

%close all
%interactiveSVDspace(U,V,numComponents,timeAxis)

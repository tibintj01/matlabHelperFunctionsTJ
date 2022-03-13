freqBandName='Alpha';
if(~exist(sprintf('singleCycle%sSVD.mat',freqBandName)))
    disp('loading single cycle matrix....')
    tic
    data=load('/nfs/turbo/lsa-ojahmed/mg49_singleCycleStats_Par/cycleMatrixMaxNormed-MG49-Alpha-zScoredCycles-100.0-ms-Window.mat');
    toc

    singleCycles=data.cycleMatrix;

    tic
    disp('computing SVD')
    [U,S,V]=svd(singleCycles(:,1:1000000),'econ');
    %[U,S,V]=svd(singleCycles(:,1:end),'econ');
    toc
    numVectorsSave=20;
    U=U(:,1:numVectorsSave);
    S=S(:,1:numVectorsSave);
    V=V(:,1:numVectorsSave);
    save(sprintf('singleCycle%sSVD.mat',freqBandName),'U','S','V');

else
    svdData=load(sprintf('singleCycle%sSVD.mat',freqBandName));
    U=svdData.U;
    S=svdData.S;
    V=svdData.V;
end

sigma=diag(S);

figure
numComponents=3;
subplot(numComponents+1,1,1)

decFs=30000/15;

plot(sigma(1:20),'r*')
xlabel('Component number')
ylabel('Singular value')
title(sprintf('%s: Strength of components across cycles',freqBandName))
timeAxis=(1:size(U,1))/decFs;
for i=1:numComponents
   subplot(numComponents+1,1,i+1)
   
   plot(timeAxis,U(:,i))
   %xlim([1 length(U(:,i))])
   xlabel('Time in cycle')
   ylabel('zscore LFP')
   title(sprintf('Scalable waveform component %d',i))
end



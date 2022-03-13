function [] = computeAndSaveSVD(freqBandName,columnDescrp,svdMatrixMatPath)
%freqBandName='Alpha';
%Description: given matrix mat file and descriptive labels, saves first 20 singular value modes describing data


if(~exist(sprintf('%s-%sSVD.mat',columnDescrp,freqBandName)))
    disp(sprintf('loading %s matrix....',columnDescrp))
    tic
    %data=load('/nfs/turbo/lsa-ojahmed/mg49_singleCycleStats_Par/cycleMatrixMaxNormed-MG49-Alpha-zScoredCycles-100.0-ms-Window.mat');
    data=load(svdMatrixMatPath);
    toc

    try
   	 A=data.A;
    catch
	A=data.cycleMatrix;
    end
    tic
    disp('computing SVD')
    maxNumCols=min(1000000,size(A,2));
    [U,S,V]=svd(A(:,1:maxNumCols),'econ');
    %[U,S,V]=svd(A(:,1:end),'econ');
    toc
    numVectorsSave=20;
    U=U(:,1:numVectorsSave);
    S=S(:,1:numVectorsSave);
    V=V(:,1:numVectorsSave);
    save(sprintf('%s-%sSVD.mat',columnDescrp,freqBandName),'U','S','V');

else
    svdData=load(sprintf('%s-%sSVD.mat',columnDescrp,freqBandName));
    U=svdData.U;
    S=svdData.S;
    V=svdData.V;
end

sigma=diag(S);

figure
numComponents=4;
subplot(numComponents+1,1,1)


plot(sigma(1:20),'r*')
xlabel('Component number')
ylabel('Singular value')
title(sprintf('%s %s: Strength of components',freqBandName,columnDescrp))

%decFs=30000/15;
%timeAxis=(1:size(U,1))/decFs;

%timeAxis=(1:size(U,1))/decFs;
for i=1:numComponents
   subplot(numComponents+1,1,i+1)
  
	%timeAxis=0.5:0.4883:119.5;
	 
   %plot(timeAxis,U(:,i))
   plot(U(:,i))
   %plot(timeAxis,log(U(:,i)))
   %xlim([1 length(U(:,i))])
   xlabel('Time in cycle')
   %xlabel('Freq (Hz)')
   ylabel('zscore LFP')
   %ylabel('ln(Power)')
   title(sprintf('Scalable waveform component %d',i))
end



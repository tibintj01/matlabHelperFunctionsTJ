function [] = concatParDec()

dataDir='C:/Users/tibin/spikeDynamicsAnalysisTibin/session3Mat';
sessNum=3;
numChan=96;

	tic
	disp('saving decimated lfp using parallel cores, 1 channel per core....')
	parfor chanNum=1:numChan
		concatenateLFPchanAcrossMatFilesV73_decChan(sessNum,dataDir,chanNum);
	end

	toc

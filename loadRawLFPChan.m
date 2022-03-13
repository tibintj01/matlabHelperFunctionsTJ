function []= loadRawLFPChan(chanNum)

	localLFPdata=matfile(sprintf('C:/Users/tibin/spikeDynamicsAnalysisTibin/session3Mat/concatChan%d/lfpOriginalConcat.mat',chanNum));

    decFactor=15;
        localLFP=localLFPdata.concatenatedLFP;
    localLFP=localLFP(1:decFactor:end);
        localLFP=getFilteredLFPLowPass(localLFP,cell.Fs/decFactor);


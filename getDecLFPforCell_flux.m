function [decLFP,decFs] = getDecLFPforCell(cellPropPath)

	cellProp=load(cellPropPath);

	chanNum=getChFromPath(cellPropPath);

        %localLFPdata=matfile(sprintf('C:/Users/tibin/spikeDynamicsAnalysisTibin/session3Mat/concatChan%d/lfpDecConcat.mat',chanNum));
        localLFPdata=matfile(sprintf('/nfs/turbo/lsa-ojahmed/spikeDynamicsAnalysisTibin/session3Mat/concatChan%d/lfpDecConcat.mat',chanNum));

	%decFs=cell.Fs/15; %make sure this dec factor holds!!!!
	decFs=30000/15; %make sure this dec factor holds!!!!

	decLFP=localLFPdata.concatenatedLFP;

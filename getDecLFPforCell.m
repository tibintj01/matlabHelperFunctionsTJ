function [decLFP,decFs] = getDecLFPforCell(cell)

	chanNum=cell.chanNum;

        localLFPdata=matfile(sprintf('C:/Users/tibin/spikeDynamicsAnalysisTibin/session3Mat/concatChan%d/lfpDecConcat.mat',chanNum));

	decFs=cell.Fs/15; %make sure this dec factor holds!!!!

	decLFP=localLFPdata.concatenatedLFP;

function [filteredLFPlowPass]=getFilteredLFPforSpikeSorting(originalLFP,Fs)
	%lowFreq=0.1/(Fs/2)
	%highFreq=300/(Fs/2) %lfpN_dec.mat is low pass filtered cutoff ~800Hz
	lowFreq=0.1;
	highFreq=300; %lfpN_dec.mat is low pass filtered cutoff ~800Hz
	filterOrder=20;
	%filterOrder=2;
	zscore=1;

	lfp=double(originalLFP);
	filteredLFPlowPass=filterLFP(lfp, Fs, lowFreq, highFreq, filterOrder, zscore);



function [filteredLFPforSpikes]=getFilteredLFPforSpikeSorting(originalLFP,Fs)
	%lowFreq=1000;
	lowFreq=300;
	%lowFreq=600;
	%lowFreq=400;
	highFreq=3000;
	%highFreq=1000;
	
	%filterOrder=20;
	filterOrder=2;
	zscore=0;

	lfp=double(originalLFP);
	filteredLFPforSpikes=filterLFP(lfp, Fs, lowFreq, highFreq, filterOrder, zscore);


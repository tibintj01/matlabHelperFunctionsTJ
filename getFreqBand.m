function [lowFreq,highFreq]=getFreqBand(freqBand)

	if(contains(freqBand,'Delta'))
		lowFreq=0.1;
		highFreq=4;
	elseif(contains(freqBand,'Alpha'))
		%lowFreq=8;
		%highFreq=12;
		lowFreq=8.2;
		highFreq=13;
	elseif(contains(freqBand,'HighGamma'))
		lowFreq=60;
		highFreq=100;
	elseif(contains(freqBand,'LowGamma'))
		lowFreq=30;
		highFreq=55;
	elseif(contains(freqBand,'LowBeta'))
		lowFreq=12.5;
		highFreq=18;
	elseif(contains(freqBand,'HighBeta'))
		lowFreq=18.5;
		highFreq=25;
	elseif(contains(freqBand,'Theta'))
		lowFreq=4.5;
		%highFreq=9.2;
		highFreq=9.5;
	end

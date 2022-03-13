function [halfPeakWidth] = getHalfPeakWidth(waveform,Fs,andPlot)
	%Tibin John, 9/6/17
	%Description: returns time difference between values at half of peak
	%Input: waveform values and sampling rate
	%Output: time between values at half of peak

	%assumes non-inverted extracellular spike around 0

	peakVal=max(waveform);
	halfPeakVal=peakVal/2;
	
	interpFact=5;
	[interpolatedWaveform,newFs]=getInterped(waveform,interpFact,Fs);

	[interpMax,interpPkIdx]=max(interpolatedWaveform);

	firstIdx=getIdxOf(halfPeakVal,interpolatedWaveform(1:interpPkIdx));

	secondIdx=getIdxOf(halfPeakVal,interpolatedWaveform(interpPkIdx:end))+interpPkIdx-1;

	halfPeakWidth=(secondIdx-firstIdx)/newFs;

	if(exist('andPlot'))
		time2=(secondIdx)/newFs
		time1=(firstIdx)/newFs
		figure
		displayWaveform(waveform,Fs)
		hold on
		
		plot([time1 time2]*1000,[halfPeakVal halfPeakVal]*1000,'k')
	end


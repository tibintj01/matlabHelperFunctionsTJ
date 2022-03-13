function [valleyToPeakTime] = getValleyPeakTimeGivenWaveform(waveform,Fs,andPlot)
	%Tibin John, 9/6/17
	%Description: Returns time between valley and peak of waveform 
	%Input: waveform values and sampling rate (plot option - pass in any variable 3rd)
	%Output: valley to peak time in seconds

	interpFact=5;
        [interpolatedWaveform,newFs]=getInterped(waveform,interpFact,Fs);

	[minVal, minIdx]=min(interpolatedWaveform);
	[maxVal,maxIdx]=max(interpolatedWaveform);

	%assumes extracellular waveforms are not inverted
	valleyToPeakTime=(maxIdx-minIdx)/newFs;

 	if(exist('andPlot'))
                time2=(maxIdx)/newFs
                time1=(minIdx)/newFs
                figure
                displayWaveform(waveform,Fs)
                hold on

                %plot([time1 time2]*1000,[maxVal maxVal]*1000,'k')
                plot([time1 time2]*1000,[0 0]*1000,'k')
        end

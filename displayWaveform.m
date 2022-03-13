function [] = displayWaveform(waveform,Fs)
	%Tibin John, 9/6/17
	%Description: displays spike waveform
	%Input: waveform values, sampling rate
	%Output: plot

	timeAxis=(1:length(waveform))/Fs;
	
	secToMs=1000;

	timeAxis=timeAxis*secToMs;

	%assumes in millivolts, converting to uV
	milliToMicro=1000;
	waveform=waveform*milliToMicro;

	plot(timeAxis,waveform,'ro')

	xlabel('Time (ms)')
	ylabel('\muV')


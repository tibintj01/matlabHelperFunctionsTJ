ch=6;
singleCycleType='alpha8_2-13Hz';
%cellLetter='a';
cellLetter='c';


plotRawLFPwithPatternIntensities(ch,singleCycleType,cellLetter)

%plotSpectrogram(decLFP,Fs,)
%getInstSpecSlopeFromLFP(decLFP,Fs)

%decLFP=decLFP(1:1e6);
%timeAxis=timeAxis(1:1e6);

%[S phi f t] = omarCoherogramLFPSpikes(timeAxis,decLFP, Fs, spikeTimes, window, winstep);

%S=S';

%figure
%imagesc(t,f,S)
%colorbar

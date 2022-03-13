function [outputStruct] = plotSpectrogram(lfp,Fs,bandNameNum,bandNameDenom,figH)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Description:
%
%
%Input:
%
%
%Output:
%
%
%Author: Tibin John, tibintj@umich.edu
%Project directory name: /nfs/turbo/lsa-ojahmed/tibin/spikeDynamicsAnalysisTibin 
%Created on 2018-06-12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%extract input variables from struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE THE SPECTROGRAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('computing output.........')
         %spectrumWindow = 1;
         spectrumWindow = 0.75;
        %spectrumWinstep = 0.5;
        %spectrumWinstep = 0.25;
        spectrumWinstep = 0.1;
        tapers = [1 1];
        %tapers = [2 3];
        %fpass = [0.5 150];
        %fpass = [0.5 30];
        fpass = [0.5 80];
        [Sgram fSgram tSgram] = omarSpectrogram(lfp, Fs, spectrumWindow, spectrumWinstep, tapers, fpass, 0);
        Sgramlog = log(Sgram);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT the spectrogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                omarPcolor(tSgram,fSgram,Sgramlog',figH);
		%omarPcolor(tSgram,fSgram,Sgramlog');
                shading flat
		cMin=prctile(Sgramlog(:),0.5);
		%cMax=prctile(Sgramlog(:),95);
		cMax=prctile(Sgramlog(:),99.5);
                caxis([cMin cMax])
                %caxis([-4 25])
                %caxis([5 25])
                %colorbar
		colormap jet;
		%colormap copper;
	 %colorbar
        %ylim([0.52 60])
        xlabel('Time (seconds)');
        ylabel('LFP Freqeuncy (Hz)');
        title('Spectrogram')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get frequency band ratio over time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(exist('bandNameNum','var') && exist('bandNameDenom','var'))
	[band1FreqLow,band1FreqHigh]=getFreqBand(bandNameNum);
	band1Indices=find(fSgram>band1FreqLow & fSgram < band1FreqHigh);
	
	[band2FreqLow,band2FreqHigh]=getFreqBand(bandNameDenom);
	band2Indices=find(fSgram>band2FreqLow & fSgram < band2FreqHigh);

	band1PowerOverTime=mean(Sgram(:,band1Indices),2);
	band2PowerOverTime=mean(Sgram(:,band2Indices),2);
	
end

outputStruct.band1PowerOverTime=band1PowerOverTime;
outputStruct.band2PowerOverTime=band2PowerOverTime;
outputStruct.bandRatio=band1PowerOverTime./band2PowerOverTime;
outputStruct.tSgram=tSgram;

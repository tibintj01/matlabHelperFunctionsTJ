function [] = plotSpectrogram(lfp,Fs,bandNameNum,bandNameDenom,figH)
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
sampPeriod=1/Fs;
lfpTimeAxis=sampPeriod:sampPeriod:(sampPeriod*length(lfp));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE THE SPECTROGRAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('computing output.........')
         %spectrumWindow = 1;
        %spectrumWinstep = 0.5;
         spectrumWindow = 0.75;
        spectrumWinstep = 0.1;
        tapers = [1 1];
        %tapers = [2 3];
        %fpass = [0.5 120];
        fpass = [0.5 200];
        [Sgram fSgram tSgram] = omarSpectrogram(lfp, Fs, spectrumWindow, spectrumWinstep, tapers, fpass, 0);
        Sgramlog = log(Sgram);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT the spectrogram and raw lfp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		subplot(3,1,1)
                omarPcolor(tSgram,fSgram,Sgramlog',figH);
		%omarPcolor(tSgram,fSgram,Sgramlog');
                shading flat
                %caxis([-4 25])
                %caxis([prctile(Sgramlog(:),10)  25])
                caxis([prctile(Sgramlog(:),10)  prctile(Sgramlog(:),95)])
                %colorbar
		colormap jet;
	 %colorbar
        %ylim([0.52 60])
        xlabel('Time (seconds)');
        ylabel('LFP Freqeuncy (Hz)');
        title('Spectrogram')


	subplot(3,1,3)
	plot(lfpTimeAxis,lfp)

        xlabel('Time (seconds)');
        ylabel('LFP (uV)');
        title('Raw LFP')

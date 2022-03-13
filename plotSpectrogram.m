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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE THE SPECTROGRAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('computing output.........')
         spectrumWindow = 1;
        spectrumWinstep = 0.5;
        tapers = [1 1];
        %tapers = [2 3];
        %fpass = [0.5 120];
        fpass = [0.5 20];
        [Sgram fSgram tSgram] = omarSpectrogram(lfp, Fs, spectrumWindow, spectrumWinstep, tapers, fpass, 0);
        Sgramlog = log(Sgram);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT the spectrogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                omarPcolor(tSgram,fSgram,Sgramlog',figH);
		%omarPcolor(tSgram,fSgram,Sgramlog');
                shading flat
                %caxis([-4 25])
                caxis([5 25])
                %colorbar
		colormap jet;
	 %colorbar
        %ylim([0.52 60])
        xlabel('Time (seconds)');
        ylabel('LFP Freqeuncy (Hz)');
        title('Spectrogram')


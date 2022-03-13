function [outputStruct] = getBandPowerOverTime(lfpTime,lfp,Fs,fBandDesired,timeWindow)
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
disp('computing spectrogram.........')
        tic
	 %spectrumWindow = 1;
	 spectrumWindow = timeWindow; %TIME WINDOW SET HERE
        spectrumWinstep = 0.5;
        tapers = [1 1];
        %tapers = [2 3];
        fpass = [0.5 120];
        %fpass = [0.5 20];
        [Sgram fSgram tSgram] = omarSpectrogram(lfp, Fs, spectrumWindow, spectrumWinstep, tapers, fpass, 0);
        %Sgramlog = log(Sgram);
	toc

	%adjust start of band time (normally 1/spectrumWinStep) so that first bin is at lfpTime(1)-1/Fs ("0") + 1/spectrumWinStep
	%specTime=tSgram+lfpTime(1)-1/Fs;
	specTime=tSgram;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get frequency band ratio over time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	band1Indices=find(fSgram>fBandDesired(1) & fSgram < fBandDesired(2));
	
	band1PowerOverTime=mean(Sgram(:,band1Indices),2);

	oneSD=std(band1PowerOverTime);
	meanPlus1SD=mean(band1PowerOverTime)+std(band1PowerOverTime);	
	meanPlus2SD=mean(band1PowerOverTime)+2*std(band1PowerOverTime);	
	meanPlus3SD=mean(band1PowerOverTime)+3*std(band1PowerOverTime);	

	%over1SDtimes=tSgram(band1PowerOverTime>=meanPlus1SD);
	%over2SDtimes=tSgram(band1PowerOverTime>=meanPlus2SD);
	%over3SDtimes=tSgram(band1PowerOverTime>=meanPlus3SD);

	over1SDtimes=specTime(band1PowerOverTime>=meanPlus1SD);
	over2SDtimes=specTime(band1PowerOverTime>=meanPlus2SD);
	over3SDtimes=specTime(band1PowerOverTime>=meanPlus3SD);

outputStruct.band1PowerOverTime=band1PowerOverTime;
outputStruct.tSgram=specTime;
outputStruct.tBinWidth=spectrumWindow;
outputStruct.tBinStep=spectrumWinstep;
outputStruct.oneSD=oneSD;

outputStruct.over1SDtimes=over1SDtimes;
outputStruct.over2SDtimes=over2SDtimes;
outputStruct.over3SDtimes=over3SDtimes;

function [C Syy Snn f] = getLfpSpikeSpectra(timeAxis, lfp, Fs, spikeTimes, tapers, fpass, avg)
% [S phi f t] = omarCoherogramLFPSpikes(timeAxis, lfp, Fs, spikeTimes, window, winstep, tapers, fpass, avg)
% DEFAULTS: tapers: [3 5], fpass = [1 250], avg = 1

if ~exist('tapers')
    %tapers = [3 5];
    tapers = [1 1];
end
if ~exist('fpass')
    %fpass = [1 250];
    %fpass = [0.5 100];
    fpass = [0.5 80];
end
if ~exist('avg')
    avg = 1;
end

%movingwin = [window winstep];

params.tapers = tapers;
params.Fs = Fs;
params.fpass = fpass;%[1 150];%[18 272];
params.pad = 0;
params.trialave = avg;

%% NOW ADJUST THE SPIKE TIMES INTO BINS
startTime = timeAxis(1);
spikeTimes = spikeTimes - startTime;


lfp=lfp(:);
spikeTimes=spikeTimes(:);

T=timeAxis(end)-timeAxis(1);

useBinnedProcess=1;
%useBinnedProcess=0;

if(useBinnedProcess)
	spikeBool=rasterize(spikeTimes,Fs,T);

	spikeBool=spikeBool(:);
	if(length(spikeBool)<length(lfp))
		spikeBool=[spikeBool; 0];
	elseif(length(spikeBool)>length(lfp))
		spikeBool=spikeBool(1:(end-1));
	end

	%size(lfp')
	%size(spikeBool')
	%[S,phi,S12,S1,S2,t,f] = cohgramcpt(lfp', spikeTimes', movingwin, params);
	%[C,~,~,Syy,Snn,f] = coherencycpb(lfp', spikeTimes',  params);
	%[C,~,~,Syy,Snn,f] = coherencycpb(lfp', spikeBool',  params);
	[C,~,~,Syy,Snn,f] = coherencycpb(lfp, spikeBool,  params);
else

	[C,~,~,Syy,Snn,f]=coherencycpt(lfp,spikeTimes,params);
end

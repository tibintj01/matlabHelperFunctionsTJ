%% LFP FILENAME
%  This .mat file contains the following variables lfp, lfpFS
lfpFilename1 = 'Psych433_Ahmed_LFP1.mat';

%% LOAD THE DATA - JUST THE ONE FILE FOR NOW
load(lfpFilename1);

%% TECHNICAL ASIDE - don't worry about this for now
nanLFPBins = find(isnan(rawLFP1));
nonNanLFPBins = find(isnan(rawLFP1) == 0);
rawLFP1(nanLFPBins) = interp1(nonNanLFPBins, rawLFP1(nonNanLFPBins), nanLFPBins, 'linear', 'extrap');


%% Fs = Sampling Rate (Fs = Frequency of sampling). Use the Fs to construct the Time Axis for the LFP
numBins = length(rawLFP1); % The number of samples (i.e. points) in our dataset
lfpTimeAxis = (1:numBins) * (1/Fs); % This converts from bins to seconds

lfpTimeAxis(1:50) % Note the lack of a semicolon here. Leaving out the semicolon will display the output of this command in the command window. Look at the command window now

% return; % "return" stops the script at this point

%% PLOT THE RAW DATA
figure;  % initialize a NEW figure
hold on; % this allows you to make multiple plots on the same axes, if you want
box off; % this is purely aesthetic - it removes the lines from the top and the right of the box
set(gca, 'fontsize', 15);
plot(lfpTimeAxis, rawLFP1);
xlabel('Time (seconds)');
ylabel('Amplitude (uV)');

%% NOW REPLOT BY SKIPPING BINS - this makes things a little faster
%  The notion of sampling rate - why does it matter?
figure;  % initialize a NEW figure
hold on; % this allows you to make multiple plots on the same axes, if you want
box off; % this is purely aesthetic - it removes the lines from the top and the right of the box
set(gca, 'fontsize', 15);
plot(lfpTimeAxis(1:5:end), rawLFP1(1:5:end));
xlabel('Time (seconds)');
ylabel('Amplitude (uV)');

%% FILTER THE DATA IN A SPECIFIC BAND - DELTA
minFiltFreq = 0.5;
maxFiltFreq = 4;
filterOrder = 2; % go with 2 most of the time
zscoreFilteredLFP = 1; % whether or not to zscore the filtered LFP
filteredLFPBand1 = filterLFP(rawLFP1, Fs, minFiltFreq, maxFiltFreq, filterOrder, zscoreFilteredLFP);

%% PLOT THE FILTERED DATA
figure;  % initialize a NEW figure
hold on; % this allows you to make multiple plots on the same axes, if you want
box off; % this is purely aesthetic - it removes the lines from the top and the right of the box
set(gca, 'fontsize', 15);
plot(lfpTimeAxis, filteredLFPBand1);
xlabel('Time (seconds)');
ylabel('Amplitude (Z)');

%% FILTER THE DATA IN ANOTHER SPECIFIC BAND - FAST BETA
minFiltFreq = 18;
maxFiltFreq = 25;
filterOrder = 2; % go with 2 most of the time
zscoreFilteredLFP = 1; % whether or not to zscore the filtered LFP
filteredLFPBand2 = filterLFP(rawLFP1, Fs, minFiltFreq, maxFiltFreq, filterOrder, zscoreFilteredLFP);


%% PLOT THE TWO FILTER BANDS ON THE SAME FIGURE USING SUBPLOT
%  What do you notice about the two bands?
figure;  % initialize a NEW figure

subplot(2,1,1);
hold on; % this allows you to make multiple plots on the same axes, if you want
box off; % this is purely aesthetic - it removes the lines from the top and the right of the box
set(gca, 'fontsize', 15);
plot(lfpTimeAxis, filteredLFPBand1, 'k');
legend('Delta Rhythms');
xlabel('Time (seconds)');
ylabel('Amplitude (Z)');
ylim([-8 8]);

subplot(2,1,2);
hold on; % this allows you to make multiple plots on the same axes, if you want
box off; % this is purely aesthetic - it removes the lines from the top and the right of the box
set(gca, 'fontsize', 15);
plot(lfpTimeAxis, filteredLFPBand2, 'r');
legend('Fast Beta Rhythms');
xlabel('Time (seconds)');
ylabel('Amplitude (Z)');
ylim([-8 8]);


%% LETS DO A POWER SPECTRUM OF THE ENTIRE DATA
%  Is this smart? What could we do better?
minFreqForSpectrum = 1;
maxFreqForSpectrum = 120;
params.Fs = Fs;
params.trialave = 1;
params.pad = 0;
params.tapers = [1 1];
[S f] = spectrumWholeSimple(rawLFP1, Fs, minFreqForSpectrum, maxFreqForSpectrum, params);


%% PLOT THE SPECTRUM of the ENTIRE DATA
figure;
hold on;
plot(f, S, 'linewidth', 1, 'color', [0.1 0.1 0.1]);
xlabel('Frequency (Hz)', 'fontsize', 10);
ylabel('Power Spectral Density (units^2 / Hz)', 'fontsize', 10);
xlim([0 50])
set(gca, 'xtick', [0:4:50]);
% set up a log y axis
set(gca, 'yscale', 'log');


%% LETS DIVIDE THE DATA INTO 2 CHUNKS - BEFORE AND AFTER SOMETHING MYSTERIOUS HAPPENS
preBins = find(lfpTimeAxis >= 0 & lfpTimeAxis < 1500); % Find the bins corresponding to the first 1500 seconds of the data
preRawLFP1 = rawLFP1(preBins);
postBins = find(lfpTimeAxis >= 3000 & lfpTimeAxis < 5000); % Find the bins corresponding to 3000-4000 seconds of the data
postRawLFP1 = rawLFP1(postBins);

% NOW LETS CALCULATE TWO POWER SPECTRA - one on each chunk
[Spre fpre] = spectrumWholeSimple(preRawLFP1, Fs, minFreqForSpectrum, maxFreqForSpectrum, params);
[Spost fpost] = spectrumWholeSimple(postRawLFP1, Fs, minFreqForSpectrum, maxFreqForSpectrum, params);

% LETS PLOT THE OVERLAID POWER SPECTRA - WHAT DOES IT SHOW?
figure;
hold on;
plot(fpre, Spre, 'linewidth', 1, 'color', [0.1 0.1 0.1]);
plot(fpost, Spost, 'linewidth', 1, 'color', [0.8 0.2 0.2]);
legend('Pre', 'Post');
xlabel('Frequency (Hz)', 'fontsize', 10);
ylabel('Power Spectral Density (units^2 / Hz)', 'fontsize', 10);
xlim([0 50])
set(gca, 'xtick', [0:4:50]);
% set up a log y axis
set(gca, 'yscale', 'log');

%% NOW CALCULATE THE SPECTROGRAM
spectrumWindow = 30;
spectrumWinstep = 30;
tapers = [1 1];
fpass = [0.5 120];
[Sgram fSgram tSgram] = omarSpectrogram(rawLFP1, Fs, spectrumWindow, spectrumWinstep, tapers, fpass, 0);
Sgramlog = log(Sgram);

%% PLOT the spectrogram
omarPcolor(tSgram,fSgram,Sgramlog');
shading flat
caxis([-4 7])
colormap jet;
colorbar
ylim([0.52 60])
xlabel('Time (seconds)');
ylabel('LFP Freqeuncy (Hz)');

%% OMAR TEMP ADDITION TO COMPARE THE COMPUTER "SPECTRUM" to TIME-LIMITED AVERAGE OF SPECTROGRAM
preSgramBins = find(tSgram >= 0 & tSgram < 1500);
sgramPreSpectrum = nanmean(Sgram(preSgramBins, :));

figure;
hold on;
plot(fpre, Spre, 'linewidth', 1, 'color', [0.1 0.1 0.1]);
plot(fSgram, sgramPreSpectrum, 'linewidth', 1, 'color', [0.2 0.8 0.2]);
legend('Pre Spectrum', 'Pre Spectrum from Spectrogram');
xlabel('Frequency (Hz)', 'fontsize', 10);
ylabel('Power Spectral Density (units^2 / Hz)', 'fontsize', 10);
xlim([0 50])
set(gca, 'xtick', [0:4:50]);
% set up a log y axis
set(gca, 'yscale', 'log');

return; % END OF MATLAB DAY 1

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MATLAB DAY 2: NOW LETS WORK WITH BOTH SPIKES AND LOCAL FIELD POTENTIAL

%% FIRST LOAD THE SPIKES
spikesFilename = 'Psych433_Ahmed_Spikes.mat';
load(spikesFilename);

%% NOW LETS LOOK AT DELTA-BAND OSCILLATIONS
minFiltFreq = 0.5;
maxFiltFreq = 4;
filterOrder = 2; % go with 2 most of the time
zscoreFilteredLFP = 1; % whether or not to zscore the filtered LFP
filtLFP = filterLFP(rawLFP1, Fs, minFiltFreq, maxFiltFreq, 2, 1);
hilbertLFP = hilbert(-filtLFP); % NOTE THE INVERSION TO GET 180=MINIMA
phaseLFP = angle(hilbertLFP);

%% Convert phase from radians to degrees
% phase ranges from -pi to +pi... make it range from 0 to 360
phaseLFP = ((phaseLFP/pi) + 1)/2 * 360;

%% TEST TO MAKE SURE IT WORKED: Let's see what the phase looks like plotted against the filtered LFP
lfpInd = find(lfpTimeAxis >= 1000.5 & lfpTimeAxis < 1001.5);
figure;
hold on;
plot(lfpTimeAxis(lfpInd), zscore(filtLFP(lfpInd)), 'r');
plot(lfpTimeAxis(lfpInd), zscore(phaseLFP(lfpInd)), 'b');

%% NOW LETS CALC THE DELTA PHASE OF EACH SPIKE
spikeBins = double(spikeBins1);
spikeTimes = spikeBins ./ spikeFs;
goodSpikeInd = find(spikeTimes < 10000); % This is just to consider the spikes in the first 10000 seconds
spikeTimes = spikeTimes(goodSpikeInd);
spikeLFPBins = round(spikeTimes * Fs); % convert from spikeBins to lfpBins
spikePhase = phaseLFP(spikeLFPBins);

%% Lets plot a histogram of these phases to see if there is a trend
phaseBinEdges = 0:20:360;
phaseBinCenters = [10:20:350];
spikeCount = histc(spikePhase, phaseBinEdges);
spikePercent = spikeCount(1:end-1) ./ sum(spikeCount) * 100;

figure;
hold on;
plot(phaseBinCenters, spikePercent, 'x-', 'linewidth', 1, 'color', [0.1 0.1 0.1]);
xlabel('Phase (degrees)');
ylabel('% of Spikes');
xlim([0 360]);

%% NOW Lets compare phase locking during the PRE and POST phases
spikePreInd = find(spikeTimes >= 0 & spikeTimes < 1500);
spikePostInd = find(spikeTimes >= 3000 & spikeTimes < 5000);
spikeCountPre = histc(spikePhase(spikePreInd), phaseBinEdges);
spikeCountPost = histc(spikePhase(spikePostInd), phaseBinEdges);
spikePercentPre = spikeCountPre(1:end-1) ./ sum(spikeCountPre) * 100;
spikePercentPost = spikeCountPost(1:end-1) ./ sum(spikeCountPost) * 100;

figure;
hold on;
plot(phaseBinCenters, spikePercentPre, 'x-', 'linewidth', 2, 'color', [0.1 0.1 0.1]);
plot(phaseBinCenters, spikePercentPost, 'o-', 'linewidth', 2, 'color', [0.8 0.2 0.2]);
xlabel('Phase (degrees)');
ylabel('% of Spikes');
xlim([0 360]);
set(gca, 'xtick', 0:60:360);
legend('Pre', 'Post');

%% Finally lets calculate the circular statistics on the pre and post data
% PRE
spikePhaseMeanPre = mod(rad2ang(circ_mean(ang2rad(spikePhase(spikePreInd)))),360);
[spikePhaseRTestPPre spikePhaseRTestZPre] = circ_rtest(ang2rad(spikePhase(spikePreInd)));
spikePhaseRPre = circ_r(ang2rad(spikePhase(spikePreInd)));
spikePhaseKappaPre   = circ_r2kappa(spikePhaseRPre);

% POST
spikePhaseMeanPost = mod(rad2ang(circ_mean(ang2rad(spikePhase(spikePostInd)))),360);
[spikePhaseRTestPPost spikePhaseRTestZPost] = circ_rtest(ang2rad(spikePhase(spikePostInd)));
spikePhaseRPost = circ_r(ang2rad(spikePhase(spikePostInd)));
spikePhaseKappaPost   = circ_r2kappa(spikePhaseRPost);

%% NOW use the command window to look at the values and compare them

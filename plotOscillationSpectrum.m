function [] = plotOscillationSpectrum(rawLFP1,Fs)

disp('computing spectrogram....')
spectrumWindow = 1;
spectrumWinstep = 1;
tapers = [1 1];
fpass = [0.5 120];
[Sgram fSgram tSgram] = omarSpectrogram(rawLFP1, Fs, spectrumWindow, spectrumWinstep, tapers, fpass, 0);
Sgramlog = log(Sgram);

%{
%% PLOT the spectrogram
omarPcolor(tSgram,fSgram,Sgramlog');
shading flat
caxis([-4 7])
colormap jet;
colorbar
ylim([0.52 60])
xlabel('Time (seconds)');
ylabel('LFP Freqeuncy (Hz)');
%}

%preSgramBins = find(tSgram >= 0 & tSgram < 1500);
preSgramBins = find(tSgram >= 0 & tSgram < max(tSgram));
sgramPreSpectrum = nanmean(Sgram(preSgramBins, :));

figure;
hold on;
%plot(fpre, Spre, 'linewidth', 1, 'color', [0.1 0.1 0.1]);
plot(fSgram, sgramPreSpectrum, 'linewidth', 1, 'color', [0.2 0.8 0.2]);
%legend('Pre Spectrum', 'Pre Spectrum from Spectrogram');
xlabel('Frequency (Hz)', 'fontsize', 10);
ylabel('Power Spectral Density (units^2 / Hz)', 'fontsize', 10);
xlim([0 50])
set(gca, 'xtick', [0:4:50]);
% set up a log y axis
set(gca, 'yscale', 'log');
 

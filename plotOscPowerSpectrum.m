function [] = plotOscillationSpectrum(rawLFP1,Fs,newMin, newMax,colorStr,useLogScale)

%% LETS DO A POWER SPECTRUM OF THE ENTIRE DATA
%  Is this smart? What could we do better?
minFreqForSpectrum = 0.5;
%maxFreqForSpectrum = 120;
%maxFreqForSpectrum = 40;
maxFreqForSpectrum = 80;
%maxFreqForSpectrum = 30;
params.Fs = Fs;
params.trialave = 1;
params.pad = 0;
params.tapers = [1 1];
tapers=[1 1];

%[S f] = spectrumWholeSimple(rawLFP1, Fs, minFreqForSpectrum, maxFreqForSpectrum, params);
[S f] = spectrumWhole(rawLFP1, Fs, minFreqForSpectrum, maxFreqForSpectrum, tapers);
	if(useLogScale)
		S=log(S);	
	end

if(exist('newMax'))
	if(length(newMax)>0)
		S=S/max(S)*newMax;
		if(length(newMin)<1)
			newMin=min(S);
		end
		S=scaledata(S,newMin,newMax);
	else
		newMax=max(S);	
	end	
end

%% PLOT THE SPECTRUM of the ENTIRE DATA
%figure;
hold on;
if(~exist('colorStr'))
plot(f, S, 'linewidth', 1, 'color', [0.1 0.1 0.1]);
else
plot(f, S, 'linewidth', 1, 'color', colorStr);

end
	if(useLogScale)
		%set(gca, 'yscale', 'log');
	end
if(~exist('newMax'))
	xlabel('Frequency (Hz)', 'fontsize', 10);
	ylabel('Power','fontsize',10)
	%ylabel('Power Spectral Density (units^2 / Hz)', 'fontsize', 10);
	%xlim([0 50])
	%set(gca, 'xtick', [0:4:50]);
	xlim([0 maxFreqForSpectrum])
	ylim([min(S) max(S)])
	%set(gca, 'xtick', [0:4:40]);
	% set up a log y axis
else
	if(length(newMin)<1)
		newMin=min(S);
	end
	ylim([newMin newMax])
	%ylabel('Power (units^2 / Hz)','fontsize',10)
	ylabel('Power')
end
	if(useLogScale)
		%set(gca, 'yscale', 'log');
	end

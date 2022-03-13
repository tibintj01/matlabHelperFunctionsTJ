function [] = histogramLog(x,edges)
	%adapted from https://www.mathworks.com/matlabcentral/answers/127988-how-to-creat-nonlinear-bin-histogram-bar-plot-with-same-bar-width

	
	%bins = linspace(floor(min(x)), ceil(max(x)), nbins);   % Define bins
xc = histcounts(x,edges);                     % Do histogram counts
figure
bins=edgesToBins(edges);
bar(log(bins), xc)                      % Plot bars against log of ‘bins’
logxts = get(gca, 'XTick')              % 'XTick' values currently logarithmic
expxts = exp(logxts);                   % Take antilog to use them as new ‘XTickLabels’
set(gca, 'XTickLabel', floor(100*expxts)/100)
	

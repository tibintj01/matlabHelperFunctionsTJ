function [colorMap] =blackGreenColormap()

	%greenColorMap = [zeros(1, 132), linspace(0, 1, 124)];
%blueColorMap = [ zeros(1, 128), linspace(0, 1, 128)];
greenColorMap = [linspace(0, 1,256)];
colorMap = [zeros(1, 256); greenColorMap; zeros(1, 256)]';

function [colorMap] =blackRedColormap()

	%greenColorMap = [zeros(1, 132), linspace(0, 1, 124)];
%blueColorMap = [ zeros(1, 128), linspace(0, 1, 128)];
redColorMap = [linspace(0, 1,256)];
colorMap = [redColorMap; zeros(1, 256); zeros(1, 256)]';

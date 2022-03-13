function [colorMap] = blackOuterColorMap()

	%greenColorMap = [zeros(1, 170), linspace(0, 1, 86)];
	greenColorMap = [zeros(1, 128), linspace(1, 0, 64),zeros(1, 64)];
	%greenColorMap = [ linspace(0, 1, 128)  , zeros(1, 256)];

%redColorMap = [linspace(1, 0, 86), zeros(1, 170)];
%redColorMap = [linspace(0, 1, 86), zeros(1, 170)];
redColorMap = [zeros(1, 64), linspace(0, 1, 64) , zeros(1, 128)];
%redColorMap = [zeros(1, 256), linspace(1, 0, 128)];

colorMap = [redColorMap; greenColorMap; zeros(1, 256)]';

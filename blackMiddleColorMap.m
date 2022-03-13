function [colorMap] = blackMiddleColorMap()

	greenColorMap = [zeros(1, 132), linspace(0, 1, 124)];
redColorMap = [linspace(1, 0, 124), zeros(1, 132)];

%colorMap = [redColorMap; greenColorMap; zeros(1, 256)]';
colorMap = [redColorMap; greenColorMap;redColorMap]'; %purple to green!
%colorMap = [greenColorMap; redColorMap; zeros(1, 256)]'; %green to red

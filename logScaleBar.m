function [] = logScaleBar(xVals,yVals)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Description:
% given xVals spaced on a log10 scale (e.g. generated by logspace) and 
% corresponding yVals, makes bar graph on log scale with equal width bars
%Input: xVals, length n
%
%
%Output: yVals, length n
%
%
%Author: Tibin John, tibintj@umich.edu
%Project directory name: /nfs/turbo/lsa-ojahmed/tibin/spikeDynamicsAnalysisTibin 
%Created on 2018-06-01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=bar(log10(xVals),yVals)

startNum=log10(xVals(1));
endNum=log10(xVals(end));
%startNum=(xVals(1));
%endNum=(xVals(end));

%set ticks in exponent space evenly
%exponents=linspace(startNum,endNum,floor(length(xVals)*2/3));
exponents=linspace(startNum,endNum,floor(length(xVals)));
%exponentsRounded=round(exponentsRounded*1000)/1000
set(gca,'Xtick',exponents)

%rename exponent values as actual values
%set(gca,'Xticklabel',num2str(10.^get(gca,'Xtick'),'%.1f'))
for i=1:length(exponents)
	tickToLabels{i}=sprintf('%.1f',10.^exponents(i));
end
tickToLabels
%set(gca,'Xticklabel',(10.^get(gca,'Xtick')))
set(gca,'Xticklabel',tickToLabels)

%xtickformat('%.1f')


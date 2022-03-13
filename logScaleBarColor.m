function [] = logScaleBar(xVals,yVals,cVals,cLims)
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
disp('creating colored bar plot in log scale....')

tic

numColorBins=64;
%threshold outside of given range!!
cVals(cVals>cLims(2))=cLims(2);
cVals(cVals<cLims(1))=cLims(1);
%cVals

colorBins=discretize(cVals,linspace(cLims(1),cLims(2),numColorBins+1));
colorBins(isnan(cVals))=-1;



colors=parula(numColorBins);

h=bar(log10(xVals),diag(yVals),'stacked');
for i=1:length(xVals)
	set(h(i),'facecolor',colors(colorBins(i),:))
end
%{
%loop through each bar to adjust color individually
for i=1:length(xVals)
	b=bar(log10(xVals(i)),yVals(i),1/length(xVals));
	%b=bar((xVals(i)),yVals(i));
	%b.Cdata(i,:)=colors(colorBins(i));
	if(colorBins(i)>0)
		b.FaceColor=colors(colorBins(i),:);
	else
		b.FaceColor=[0 0 0];
	end
	%b.Cdata(i,:)=colors(colorBins(i));
	hold on
end
%colorbar
caxis(cLims)
%}


startNum=log10(xVals(1));
endNum=log10(xVals(end));
%startNum=(xVals(1));
%endNum=(xVals(end));

%set ticks in exponent space evenly
exponents=linspace(startNum,endNum,floor(length(xVals)*2/3));
%exponents=linspace(startNum,endNum,floor(length(xVals)));
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
xlim([-Inf Inf])
%{
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
%}

toc


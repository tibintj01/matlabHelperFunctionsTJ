
close all; clear all; clc
setTightSubplots_Low

xLog=(0:0.001:1);


base=sqrt(2)
base=1.25
base=1.3333;
base=1.4;
%yLog=getLogTime(1-xLog,base);
yLog=getLogTime(1-xLog,base);
yLog=scaledata(yLog,0,1);

figure

subplot(2,1,1)
plot(xLog,yLog,'k-','LineWidth',3)
box off
%ylim([0 1])
%xlim([0 1])

xlabel('x (a.u.)')
ylabel(sprintf('y=log_{%.2f}(1-x) (a.u.)',base))

title(sprintf('y=log_{%.2f}(1-x) (a.u.)',base))


ylim([0 1])
xlim([0 1])


box off
hold on

%plot even vertical bins to predicted bin widths
numBins=7;
vertEdges=linspace(0,1,numBins+1);
vertBinCenters=edgesToBins(vertEdges);

xEdgeVals=NaN(numBins+1,1);
xBinCenters=NaN(numBins,1);


for ei=1:numBins+1
    currEdgeVal=vertEdges(ei);
    
    
    [~,currValIdx]=min(abs(currEdgeVal-yLog));
    
    currValX=xLog(currValIdx);
    
    
    plot( [1 currValX], [ currEdgeVal currEdgeVal ],'k--','LineWidth',3)
    hold on
    plot([currValX currValX],[0 currEdgeVal],'r--','LineWidth',3)
    
    xEdgeVals(ei)=currValX;
    if(ei>1)
        xBinCenters(ei-1)=(xEdgeVals(ei-1)+xEdgeVals(ei))/2;
    end
end

xBinCenters=xBinCenters(1:end);
xEdgeVals=xEdgeVals(1:end);

plot(xBinCenters,zeros(size(xBinCenters)),'r.','MarkerSize',60)
subplot(2,1,2)
%plot(xBinCenters,diff(xEdgeVals),'k-o','LineWidth',3)


binWidths=abs(diff(xEdgeVals));
plot(xBinCenters,binWidths,'r.','MarkerSize',60)
hold on
[m,b,R]=getLinearFit(xBinCenters,binWidths,0,1)
xlim([0 1])
ylim([0 0.4])
ylabel('x bin width')
xlabel('x (a.u.)')

title({sprintf('m=%.3f',m),sprintf('e^{-m}=%.3f',exp(-m))})
box off
setFigFontTo(16)

%saveEPS('recoverLogBaseFromLinearVarIncrease1.eps')
saveas(gcf,sprintf('recoverLogBaseFromLinearVarIncrease%3f.png',base))



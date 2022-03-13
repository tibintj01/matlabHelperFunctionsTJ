
close all;clear all; clc
setTightSubplots_Low

xLog=0:0.001:1;



base=1.25

base=sqrt(2)
base=2
base=1.3333;

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
ylabel(sprintf('y=log_{%.2f}(x) (a.u.)',base))

title(sprintf('y=log_{%.2f}(x) (a.u.)',base))


ylim([0 1])
xlim([0 1])


box off
hold on

%plot even vertical bins to predicted bin widths
numBins=7;
%numBins=20;
horzEdges=linspace(0,1,numBins+1);
xBinCenters=edgesToBins(horzEdges);

yEdgeVals=NaN(numBins+1,1);


for ei=1:numBins+1
    currEdgeVal=horzEdges(ei);
    
    
    [~,currValIdx]=min(abs(currEdgeVal-xLog));
    
    currValY=yLog(currValIdx);
    
    
    plot( [ currEdgeVal currEdgeVal ],[0 currValY], 'k--','LineWidth',3)
    hold on
    plot([1 currEdgeVal],[currValY currValY],'r--','LineWidth',3)
    
    yEdgeVals(ei)=currValY;
    if(ei>1)
        yBinCenters(ei-1)=(yEdgeVals(ei-1)+yEdgeVals(ei))/2;
    end
end

yBinCenters=yBinCenters(2:end);
yEdgeVals=yEdgeVals(2:end);
xBinCenters=xBinCenters(2:end);

plot(ones(size(yBinCenters)),yBinCenters,'r.','MarkerSize',60)
subplot(2,1,2)
%plot(xBinCenters,diff(xEdgeVals),'k-o','LineWidth',3)


%{
plot(yBinCenters,diff(yEdgeVals),'r.','MarkerSize',60)
hold on
[m,b,R]=getLinearFit(yBinCenters,diff(yEdgeVals),0,1)
%}

plot(xBinCenters,diff(yEdgeVals),'r.','MarkerSize',60)
hold on
[m,b,R]=getLinearFit(xBinCenters,diff(yEdgeVals),0,1)
xlim([0 1])
%ylim([0 1])
ylabel('y bin width')
xlabel('y (a.u.)')

title({sprintf('m=%.3f',m),sprintf('e^m=%.3f',exp(m))})

setFigFontTo(16)

%saveEPS('recoverLogBaseFromLinearVarIncrease1.eps')
saveas(gcf,sprintf('recoverLogBaseFromLinearVarDecrease%3f.png',base))



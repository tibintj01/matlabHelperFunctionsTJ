
%close all; clear all; clc

xLog=0:0.001:1;


%base=sqrt(2);
%base=3;





bases=[1.4 2 3];
bases=[1.2 1.4 2 3];

%bases=[1.2 1.4 1.7 3];
numBases=length(bases);

figure

for bi=1:numBases
    base=bases(bi)
    yLog=getLogTime(1-xLog,base);
    yLog=scaledata(yLog,0,1);

    subplot(1,numBases,bi)
plot(xLog,yLog,'k-','LineWidth',3)
box off
%ylim([0 1])
%xlim([0 1])

xlabel('x (a.u.)')
ylabel('y=log(1-x) (a.u.)')

title(sprintf('log_{%.2f}(1-x)',base))

ylim([0 1])
xlim([0 1])



box off
hold on

%numBins=9;
numBins=7;
xEdges=linspace(0,1,numBins+1);
%xEdges=1-0.1*1.4.^[0 1 2 3 4 5 6 7];

xBinCenters=edgesToBins(xEdges);

yEdgeVals=NaN(numBins+1,1);
yBinCenters=NaN(numBins,1);


for ei=1:numBins+1
    currEdgeVal=xEdges(ei);
    
    
    [~,currValIdx]=min(abs(currEdgeVal-xLog));
    
    currValY=yLog(currValIdx);
    
    
    plot(  [ currEdgeVal currEdgeVal ],[0 currValY],'k--','LineWidth',3)
    hold on
    plot([1 currEdgeVal],[currValY currValY],'r--','LineWidth',3)
    
    yEdgeVals(ei)=currValY;
    if(ei>1)
        yBinCenters(ei-1)=(yEdgeVals(ei-1)+yEdgeVals(ei))/2;
    end
end

%xBinCenters=xBinCenters(1:end);
%xEdgeVals=xEdgeVals(1:end);


axis tight
box off

end

setFigFontTo(18)

saveEPS('logBaseInterpretation.eps')



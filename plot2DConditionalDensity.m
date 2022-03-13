function [] = plot2DConditionalDensity(data3cols)
%plots 2D probability density of specified subset relative to unspecified
%data
setTightSubplots_SpaceTime
%phaseEnds=(1:18)*20

numRanges=18;
phaseEnds=(numRanges:-1:3)*20
%phaseEnds=linspace(360,80,18)
numRows=2;

phaseStarts=phaseEnds-median(abs(diff(phaseEnds)));
    fH=figure;

for i=1:length(phaseStarts)
    zRangeOfInterest=[phaseStarts(i) phaseEnds(i)];
    minDataCount=20;
    minDataCount=2;
    minDataCount=5;
    minDataCount=10;
        minDataCount=6;
    xVals=data3cols(:,1);
    yVals=data3cols(:,2);
    zVals=data3cols(:,3);

    zIdxesOfInterest=zVals>=zRangeOfInterest(1) & zVals<=zRangeOfInterest(2);
    zIdxesOfInterest=logical(zIdxesOfInterest);

    manualXMax=1;
    manualYMax=1;

    histBinWidth=0.025;
    histBinWidth=0.01;
    histBinWidth=0.02;
    %histBinWidth=0.04;
    %xEdges=min(xVals):histBinWidth:max(xVals);
    %yEdges=min(yVals):histBinWidth:max(yVals);
    xVals(xVals>manualXMax)=NaN;
    yVals(yVals>manualYMax)=NaN;

    xVals(isnan(zVals))=NaN;
    yVals(isnan(zVals))=NaN;

    xEdges=min(xVals):histBinWidth:manualXMax;
    yEdges=min(yVals):histBinWidth:manualYMax;

    [pXYn,pXYc]=hist3([xVals(:),yVals(:)],{xEdges,yEdges});

    pXYn(pXYn<minDataCount)=0;

    xAndZvals=xVals;
    yAndZvals=yVals;

    xAndZvals(~zIdxesOfInterest)=NaN;
    yAndZvals(~zIdxesOfInterest)=NaN;

    [pXYZn,pXYZc]=hist3([xAndZvals(:),yAndZvals(:)],{xEdges,yEdges});

    pXYZn(pXYZn<minDataCount)=0;

    pXY=pXYn/sum(pXYn(:));
    pXYZ=pXYZn/sum(pXYZn(:));

   


    pXY=nanGaussSmoothPhaseProb(pXY);
     pXYgivenZ=pXYZ./pXY;
    pXYgivenZ=nanGaussSmoothPhaseProb(pXYgivenZ);

pXYgivenZ(pXYn<minDataCount)=NaN;
    subplot(numRows,numRanges/numRows,i)
    %omarPcolor(pXYZc{1},pXYZc{2},pXYZ',fH)
        omarPcolor(pXYZc{1},pXYZc{2},pXYgivenZ',fH)
    colormap(gca,jet)

    title(sprintf('p(X,T|%d<P<%d)',round(zRangeOfInterest(1)),round(zRangeOfInterest(2))))
    daspect([1 1 1])
    %caxis([0 0.01])
    caxis([0 4])
    drawnow
    %xlabel('X (field frac)')
    %ylabel('T (field frac)')
    
    
end
maxFig
%setFigFontTo(18)
%{
subplot(3,1,1)
omarPcolor(pXYc{1},pXYc{2},pXY',fH)
colormap(gca,jet)
cb=colorbar
ylabel(cb,'probability density')
title('p(X,T)')
%pXYgivenZ(pXYgivenZ==0)=NaN;
daspect([1 1 1])

subplot(3,1,2)
omarPcolor(pXYZc{1},pXYZc{2},pXYZ',fH)
colormap(gca,jet)
cb=colorbar
ylabel(cb,'p(X,T,P)')
title(sprintf('p(X,T,%d<P<%d)',zRangeOfInterest(1),zRangeOfInterest(2)))
daspect([1 1 1])

subplot(3,1,3)
omarPcolor(pXYZc{1},pXYZc{2},pXYgivenZ',fH)
colormap(gca,jet)
cb=colorbar
ylabel(cb,'liklihood')
title(sprintf('p(X,T | %d<P<%d)',zRangeOfInterest(1),zRangeOfInterest(2)))
daspect([1 1 1])
%}

end


function [rhoSpeedChange,slopeSpeedChange,pSpeedChange, rhoDist,slopeDist,pDist] = getSpeedChangeVsPositionSaturatingFit(data)
%fits parabola to speed vs position to get ds/dx

showPlots=0;

trackLength=120;

maxSpeedChange=5; %cm/s per cm
maxSpeedChange=3.5; %cm/s per cm
%maxSpeedChange=10; %cm/s per cm
startSaturatingSpeed=2.9;
maxReasonableSpeedChange=2.5;
maxReasonableSpeedChange=2;

minMinorAxis=45; %start at least by 60-45=15cm

speedVsPos=data.spikePerCycleInfo.speedPerCycle;
positionsCm=data.spikePerCycleInfo.normPosPerCycle*trackLength;

animalID=getFileNameFromPath(data.spikePerCycleInfo.filePath);
fieldStart=data.manualFieldStartCm;
fieldEnd=data.manualFieldEndCm;

fitEllipsedWithFixedBounds=1;
minSpeed=5; %cm/s
%minSpeed=0; %cm/s
%minSpeed=10; %cm/s
%minSpeed=15; %cm/s
%minSpeed=2; %cm/s
%minSpeed=0; %cm/s


cleanPointsFirst=1;

%minSpeed=7; %cm/s

sufficientSpeed=speedVsPos>minSpeed;

notNaNSpeed=~isnan(speedVsPos);
notNaNPositions=~isnan(positionsCm);

goodIdxes=sufficientSpeed & notNaNSpeed & notNaNPositions;

speedVsPos=speedVsPos(goodIdxes);
positionsCm=positionsCm(goodIdxes);

numGoodDataPts=length(positionsCm);

speedEdges=min(speedVsPos):max(speedVsPos);
positionEdges=min(positionsCm):max(positionsCm);

speedBins=edgesToBins(speedEdges);
posBins=edgesToBins(positionEdges);


[speedPosJoint,~,~]=histcounts2(positionsCm,speedVsPos,positionEdges,speedEdges);

if(cleanPointsFirst)
   

    smoothedSpeedPosJoint=imgaussfilt(speedPosJoint,2);
    %smoothedSpeedPosJoint=imgaussfilt(speedPosJoint,1);
    %smoothedSpeedPosJoint=imgaussfilt(smoothedSpeedPosJoint,5);
    imThresh=graythresh(smoothedSpeedPosJoint(:))*1;
    cleansedMask=smoothedSpeedPosJoint>imThresh;
    cleansedMask=imclose(cleansedMask,strel('disk',10));
    if(showPlots)
         figure; 
        subplot(3,1,1);
        imagesc(flipud(speedPosJoint'))
        colorbar

        subplot(3,1,2);
        imagesc(flipud(smoothedSpeedPosJoint'))
        colorbar

        subplot(3,1,3)

        imagesc(flipud(cleansedMask'))
        colorbar
    end

%{
numPosBins=length(posBins);
comSpeedPerPos=NaN(numPosBins,1);
for bi=1:numPosBins
    comSpeedPerPos(bi)=getCOM(smoothedSpeedPosJoint(bi,:));
end
plot(posBins,comSpeedPerPos,'ko')
hold on
plot(posBins,comSpeedPerPos,'k-')
%}


%filter by cleaned mask
keepIDs=ones(numGoodDataPts,1);
for i=1:numGoodDataPts
    currPos=positionsCm(i);
    currSpeed=speedVsPos(i);
    
    %find closest 2D bin and check if it passes mask
    [~,closestPosBinID]=min(abs(posBins-currPos));
    [~,closestSpeedBinID]=min(abs(speedBins-currSpeed));
    
    if(~cleansedMask(closestPosBinID,closestSpeedBinID))
        keepIDs(i)=0;
    end
end

cleanIDs=logical(keepIDs);

cleanPositionsCm=positionsCm(cleanIDs);
cleanSpeedVsPos=speedVsPos(cleanIDs);
end

if(showPlots)
    fH=figure; 
    hold on
    %subplot(4,1,4)
    %plot(posBins,speedFit,'b-','LineWidth',5)
    plot(positionsCm,speedVsPos,'k.')
    hold on
end

fitParabola=0;
if(fitParabola)
    %fit 2nd degree polynomial
    [pFit,S]=polyfit(cleanPositionsCm,cleanSpeedVsPos,2);
    speedFit=polyval(pFit,posBins);
else
    %ellipse_t = fit_ellipse(cleanPositionsCm,cleanSpeedVsPos,fH)
    if(cleanPointsFirst)
        ellipse_t = fit_ellipse(cleanPositionsCm,cleanSpeedVsPos);
    else
         ellipse_t = fit_ellipse(positionsCm,speedVsPos);
    end
    
    if(fitEllipsedWithFixedBounds)
     ellipseCenter=[trackLength/2; 0];
     centerPositionIDs=cleanPositionsCm>55 & cleanPositionsCm<65;
     fixedAxisHeight=nanmean(cleanSpeedVsPos(centerPositionIDs));
     [bestW,bestH]=fitHalfEllipseAxisLengths(ellipseCenter,cleanPositionsCm,cleanSpeedVsPos,fixedAxisHeight);
       bestW=max(minMinorAxis,bestW);
    end
end
if(showPlots)
    hold on
    xlim([0 max(positionEdges)])
    ylim([minSpeed max(speedEdges)])
    if(cleanPointsFirst)
        plot(cleanPositionsCm,cleanSpeedVsPos,'bo','MarkerSize',3)
    end
    hold on
end
%plotEllipse(ellipse_t)

%ellipse_t.phi=0; %assume track symmetry
%ellipse_t.X0=trackLength/2; %assume track symmetry

if(~(fitEllipsedWithFixedBounds))
    minMajorAxis=61;
    minMajorAxis=0;
    %minMajorAxis=62;
    %minMajorAxis=63;
    ellipse_t.a=max(ellipse_t.a,minMajorAxis);
    verticalShiftConstraint=10;
    verticalStretch=20;
    verticalShiftConstraint=15;
    verticalStretch=30;
    verticalShiftConstraint=30;
    verticalStretch=60;
    verticalStretch=0;
    verticalShiftConstraint=0;
    ellipse_t.b=ellipse_t.b+ellipse_t.Y0-verticalShiftConstraint+verticalStretch;
    ellipse_t.Y0=-verticalShiftConstraint;
else
    ellipse_t.phi=0;
    ellipse_t.a=bestW;
    ellipse_t.b=fixedAxisHeight;
    ellipse_t.X0=ellipseCenter(1);
    ellipse_t.Y0=ellipseCenter(2);
end

[ellipseX,ellipseSpeeds]=plotEllipse(ellipse_t);

if(showPlots)
    ylim([0 Inf])
    daspect([1 1 1])
end

%figure;
positiveSpeedIDs=ellipseSpeeds>=0;

ellipseSpeeds=ellipseSpeeds(positiveSpeedIDs);
ellipseX=ellipseX(positiveSpeedIDs);

[~,sortedIDs]=sort(ellipseX);
ellipseSpeeds=ellipseSpeeds(sortedIDs);
ellipseX=ellipseX(sortedIDs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%add saturating extrapolation up to minSpeed of fitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
goodFitIdxes=ellipseSpeeds>=minSpeed;
ellipseX=ellipseX(goodFitIdxes);
ellipseSpeeds=ellipseSpeeds(goodFitIdxes);


resampledX=linspace(min(ellipseX),max(ellipseX),length(ellipseX));
resampledSpeeds=interp1(ellipseX,ellipseSpeeds,resampledX);

delX=median(diff(resampledX)); %cm

speedSpaceDiff=diff(resampledSpeeds)/delX;
speedSpaceDiffX=resampledX(1:(end-1));
sampPointsX=linspace(0,trackLength,10000);
reasonableSpeedIdxes=abs(speedSpaceDiff)<=startSaturatingSpeed;
speedSpaceDiff=speedSpaceDiff(reasonableSpeedIdxes);
speedSpaceDiffX=speedSpaceDiffX(reasonableSpeedIdxes);

speedSpaceDiffX=[0;speedSpaceDiffX(:);trackLength];
speedSpaceDiff=[maxSpeedChange;speedSpaceDiff(:); -maxSpeedChange];

resampSpeedSpaceDiff=interp1(speedSpaceDiffX,speedSpaceDiff,sampPointsX,'pchip','extrap');
delXSaturating=median(diff(sampPointsX));

speedChangeDists=sampPointsX;
speedChanges=resampSpeedSpaceDiff;
maxReasonableSpeedCorrespondingDist=interp1(speedChanges,speedChangeDists, maxReasonableSpeedChange);

convertDistanceToSpeedChange=0;
%%%%%%%%%%%%%%%%%%%%%%%%%
%make plots
%%%%%%%%%%%%%%%%%%%%%%%%%
if(showPlots)
    
    figure;
    subplot(2,2,1)
    plot(resampledX,resampledSpeeds,'r','LineWidth',5)
    hold on
    plot(positionsCm,speedVsPos,'k.')
    xlabel('Distance along track (cm)')
    ylabel('Speed (cm/s)')
    title({'Speed vs position (semi-ellipse fit)'})
    ylim([0 max(resampledSpeeds)*1.15])
    xlim([0 trackLength])

    plot([fieldStart fieldStart],[0 max(resampledSpeeds)*1.15],'m--')
    plot([fieldEnd fieldEnd],[0 max(resampledSpeeds)*1.15],'m--')

    subplot(2,2,3)
    


    plot(sampPointsX,resampSpeedSpaceDiff,'b','LineWidth',5)



    %ylim([-5 5])
    yLims=[-3 3];
    ylim(yLims)
    xlabel('Distance along track (cm)')
    ylabel('Speed change with distance (cm/s per cm)')
    xlim([0 trackLength])
    title('Speed change vs position (semi-ellipse fit)')
    hold on

    plot([fieldStart fieldStart],yLims,'m--')
    plot([fieldEnd fieldEnd],yLims,'m--')


    subplot(2,2,2)



end
    plotPhaseVsDistance
    
if(showPlots)
    xlim([fieldStart fieldEnd])

    title({'Phase vs distance',sprintf('circlinear: rho = %.4f, slope= %.2f, p=%.2f',circLinRho,circLinSlope,circLinP)})
    subplot(2,2,4)
end



rhoDist= circLinRhoDist;
pDist=circLinPDist;
slopeDist=circLinSlopeDist;


convertDistanceToSpeedChange=1;

speedChangeBounds=interp1(speedChangeDists,speedChanges,[fieldStart fieldEnd]);
realisticSpeedChangeIdxes=abs(speedChanges)<=maxReasonableSpeedChange;
speedChanges=speedChanges(realisticSpeedChangeIdxes);
speedChangeDists=speedChangeDists(realisticSpeedChangeIdxes);



%speedChangeBounds=[min(speedChanges) max(speedChanges)];
plotPhaseVsDistance

rhoSpeedChange= circLinRhoSpeed;
pSpeedChange=circLinPSpeed;
slopeSpeedChange=circLinSlopeSpeed;

%absoluteLowerBound=max(min(speedChangeBounds),prctile(allSpeedChanges,10));
%absoluteUpperBound=min(max(speedChangeBounds),prctile(allSpeedChanges,90));


%make distance limits comparable by backward interpolation


%minSpeedChange=max(min(speedChangeBounds),absoluteLowerBound);
%maxSpeedChange=min(max(speedChangeBounds),absoluteUpperBound);

%minSpeedChange=min(speedChangeBounds);
%maxSpeedChange=max(speedChangeBounds);
%{
if(isnan(minSpeedChange))
    minSpeedChange=min(currSpeedChanges);
end

if(isnan(maxSpeedChange))
    maxSpeedChange=max(currSpeedChanges);
end

if(minSpeedChange>=maxSpeedChange)
    maxSpeedChange=3;
end
%}

if(showPlots)
    try
        title({'Phase vs speed change',sprintf('circlinear: rho = %.4f, slope= %.2f, p=%.2f',circLinRhoSpeed,circLinSlopeSpeed,circLinPSpeed)})
    catch
        title({'Phase vs speed change',sprintf('circlinear: rho = %.4f, slope= %.2f, p=%.2f',circLinRhoDist,circLinSlopeDist,circLinPDist)})
    end
    xlim([-2 2])
    hold on
    %plot(-[speedChangeBounds(1) speedChangeBounds(1)],[0 360],'m--')
    %plot(-[speedChangeBounds(2) speedChangeBounds(2)],[0 360],'m--')

    xlim(-speedChangeBounds)



    %distBounds=interp1(speedChanges,speedChangeDists,[absoluteUpperBound absoluteLowerBound ],'linear','extrap');
    subplot(2,2,2)
    %xlim(distBounds)
    xlim([0 120])
    plot([fieldStart fieldStart],[0 360],'m--')
    plot([fieldEnd fieldEnd],[0 360],'m--')
    xlim([fieldStart fieldEnd])

    setFigFontTo(18)
    uberTitle(removeUnderscores(animalID))
end




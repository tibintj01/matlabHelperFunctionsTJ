function [ellipse_t] = getSpeedChangeVsPosition(data)
%fits parabola to speed vs position to get ds/dx

trackLength=120;
speedVsPos=data.spikePerCycleInfo.speedPerCycle;
positionsCm=data.spikePerCycleInfo.normPosPerCycle*trackLength;

animalID=getFileNameFromPath(data.spikePerCycleInfo.filePath);
fieldStart=data.manualFieldStartCm;
fieldEnd=data.manualFieldEndCm;

minSpeed=5; %cm/s
minSpeed=2; %cm/s
%minSpeed=0; %cm/s


cleanPointsFirst=0;

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
    figure; 
    subplot(3,1,1);
    imagesc(flipud(speedPosJoint'))
    colorbar

    smoothedSpeedPosJoint=imgaussfilt(speedPosJoint,5);
    %smoothedSpeedPosJoint=imgaussfilt(speedPosJoint,1);
    %smoothedSpeedPosJoint=imgaussfilt(smoothedSpeedPosJoint,5);

    subplot(3,1,2);
    imagesc(flipud(smoothedSpeedPosJoint'))
    colorbar

    subplot(3,1,3)
    imThresh=graythresh(smoothedSpeedPosJoint(:))*1;
    cleansedMask=smoothedSpeedPosJoint>imThresh;
    cleansedMask=imclose(cleansedMask,strel('disk',15));
    imagesc(flipud(cleansedMask'))
    colorbar

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


fH=figure; 
hold on
%subplot(4,1,4)
%plot(posBins,speedFit,'b-','LineWidth',5)
plot(positionsCm,speedVsPos,'k.')
hold on
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
     %ellipseCenter=[trackLength/2; 0];
     %[bestW,bestH]=fitHalfEllipseAxisLengths(ellipseCenter,cleanPositionsCm,cleanSpeedVsPos);
end

hold on
xlim([0 max(positionEdges)])
ylim([minSpeed max(speedEdges)])
if(cleanPointsFirst)
    plot(cleanPositionsCm,cleanSpeedVsPos,'bo','MarkerSize',3)
end
hold on
%plotEllipse(ellipse_t)

ellipse_t.phi=0; %assume track symmetry
ellipse_t.X0=trackLength/2; %assume track symmetry

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

%verticalStretch=0;
%verticalShiftConstraint=3;


ellipse_t.b=ellipse_t.b+ellipse_t.Y0-verticalShiftConstraint+verticalStretch;
ellipse_t.Y0=-verticalShiftConstraint;

[ellipseX,ellipseSpeeds]=plotEllipse(ellipse_t);

ylim([0 Inf])

daspect([1 1 1])

for li=2:2
    %fH.Children.Children(li).LineWidth=5;
end

%figure;
positiveSpeedIDs=ellipseSpeeds>=0;

ellipseSpeeds=ellipseSpeeds(positiveSpeedIDs);
ellipseX=ellipseX(positiveSpeedIDs);

[~,sortedIDs]=sort(ellipseX);
ellipseSpeeds=ellipseSpeeds(sortedIDs);
ellipseX=ellipseX(sortedIDs);

resampledX=linspace(min(ellipseX),max(ellipseX),length(ellipseX));

resampledSpeeds=interp1(ellipseX,ellipseSpeeds,resampledX);

delX=median(diff(resampledX)); %cm

%%%%%%%%%%%%%%%%%%%%%%%%%
%make plots
%%%%%%%%%%%%%%%%%%%%%%%%%
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
speedSpaceDiff=diff(resampledSpeeds)/delX;
speedSpaceDiffX=resampledX(1:(end-1));
sampPointsX=linspace(0,trackLength,100);
resampSpeedSpaceDiff=interp1(speedSpaceDiffX,speedSpaceDiff,sampPointsX,'linear','extrap');
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



convertDistanceToSpeedChange=0;

plotPhaseVsDistance

xlim([fieldStart fieldEnd])

title('Phase vs distance')


subplot(2,2,4)

convertDistanceToSpeedChange=1;
speedChangeDists=resampledX(1:(end-1));
speedChanges=diff(resampledSpeeds)/delX;

speedChangeBounds=interp1(speedChangeDists,speedChanges,[fieldStart fieldEnd],'linear','extrap');
plotPhaseVsDistance

%absoluteLowerBound=max(min(speedChangeBounds),prctile(allSpeedChanges,10));
%absoluteUpperBound=min(max(speedChangeBounds),prctile(allSpeedChanges,90));

absoluteLowerBound=prctile(allSpeedChanges,10);
absoluteUpperBound=prctile(allSpeedChanges,90);
if(absoluteUpperBound>3)
    absoluteUpperBound=3;
end
if(absoluteLowerBound<-3)
    absoluteLowerBound=-3;
end

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

try
xlim([absoluteLowerBound absoluteUpperBound])
catch
    xlim([-5 5])
end
title('Phase vs speed change')
xlim([-2 2])
hold on
%plot(-[speedChangeBounds(1) speedChangeBounds(1)],[0 360],'m--')
%plot(-[speedChangeBounds(2) speedChangeBounds(2)],[0 360],'m--')
try
xlim(-speedChangeBounds)
end

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


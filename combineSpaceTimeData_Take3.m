%close all; clear all; clc


plotIndividualPhasePrecession=0;

for dummyRunI=2:2
    if(dummyRunI==1)
        dummyRun=0;
    else
        dummyRun=1;
    end
dataDir='/Users/tibinjohn/thetaSeq/code/spaceTimeProcessedDir';
%setTightSubplots_SpaceTime


combineDirections=1;
%combineDirections=0;

normalizeEachFieldPhaseRange=1;
normalizeEachFieldPhaseRange=0;

minPhaseDisp=120;
maxPhaseDisp=250;


for runDir=1:2
    if(runDir==1)
        trackDir=1;
    else
        trackDir=2;
    end

filePaths=getFilePathsRegex(dataDir,sprintf('*dir%d*mat',trackDir));

numSpaceBins=50;
%numSpaceBins=25;
%numSpaceBins=15;
numTimeBins=15;
%numTimeBins=10;

trackLength=120; %cm
%minSpaceCorr=-1;
%maxSpaceCorr=-0.5;

minSpaceCorr=0.35;
minSpaceCorr=0.3;
%minSpaceCorr=0.5;
%minSpaceCorr=0.3;
%minSpaceCorr=0.25;
%minSpaceCorr=0.5;
maxSpaceCorr=1;

%minSpaceCorr=-0.3
%minSpaceCorr=0.6
numFilesToUse=length(filePaths);
%numFilesToUse=100;
totalFieldResponseAcrossCells=NaN(numSpaceBins,numTimeBins,numFilesToUse);


%for i=1:length(filePaths)
%for i=1:9
skipCount=0;
if(combineDirections==1)
    if(runDir==1)
        allFieldSpaceTimePhaseTriplets=[];
        numCellsUsed=0;
    end
else
    allFieldSpaceTimePhaseTriplets=[];
    numCellsUsed=0;
end


 plotIndividualCells=0;
 if(plotIndividualCells)
      fT=figure;
 end
 minNumSpaceBinsIndiv=10;
  minNumTimeBinsIndiv=10;
  minSpeed=0.05;
  maxSpeed=Inf;

corrPThresh=0.05;
corrPThresh=0.01;
minNumSpikes=100;
numDensityBins=30;
allCellExcessJoint=zeros(numDensityBins,numDensityBins);
for i=1:numFilesToUse
    
    currFilePath=filePaths{i};
    currData=load(currFilePath);
    
    %spaceVals=currData.spaceTimePhaseInfo.spaceBins;
    %timeVals=currData.spaceTimePhaseInfo.timeBins;
    %if(isempty(spaceVals) || sum(isnan(spaceVals)) >0)
    %    continue
    %end

    linearPos=currData.spaceTimePhaseInfo.spaceTimePhasePerCycle(:,1);
    phasesInRadians=currData.spaceTimePhaseInfo.spaceTimePhasePerCycle(:,3)*pi/180;


    %if(isempty(linearPos) || isempty(phasesInRadians))
    if(length(linearPos)< minNumSpikes || length(phasesInRadians) < minNumSpikes)
	skipCount=skipCount+1;
        continue
    end

    [ rhoK,pK,sK,bK ] = kempter_lincirc( linearPos,phasesInRadians);

    showInfoPlots=0;
    [IposPhase,excessJointSpikeDensity,Hpos,Hphase,rhoSpear,pSpear,dispXlims,dispYlims]=getDiscreteMutualInfo(linearPos,phasesInRadians,numDensityBins,numDensityBins,'position','phase (radians)',showInfoPlots);

    %if(pK>corrPThresh || isnan(pK))
    if(pSpear>corrPThresh || isnan(pK) || rhoSpear>0 || isnan(rhoSpear))
	skipCount=skipCount+1;
        continue  
    end
    absSpaceTimePhaseTriplets=currData.spaceTimePhaseInfo.spaceTimePhasePerCycle;
    numDataCycles=size(absSpaceTimePhaseTriplets,1);
    


    fieldRelativeSpace=scaledata(absSpaceTimePhaseTriplets(:,1),0,1);
    
    %if(combineDirections==1 && runDir==1)
    %if(runDir==1)
    if(runDir==2)
        fieldRelativeSpace=1-fieldRelativeSpace;
    end
    fieldRelativeTime=scaledata(absSpaceTimePhaseTriplets(:,2),0,1);
    
    if(normalizeEachFieldPhaseRange)
       absSpaceTimePhaseTriplets(:,3)=scaledata(absSpaceTimePhaseTriplets(:,3),0,360);
    end
    
    % spaceBins=scaledata(currData.spaceTimePhaseInfo.spaceBins,0,1);
    %if(combineDirections==1 && runDir==1)
    %if(runDir==1)
    %if(runDir==2)
    %    spaceBins=1-spaceBins;
    %end
   %timeBins=scaledata(currData.spaceTimePhaseInfo.timeBins,0,1);
   
         %estSpeed(numCellsUsed+1)=range(currData.spaceTimePhaseInfo.spaceBins)/range(currData.spaceTimePhaseInfo.timeBins);

   
         
    %currFieldSpaceTimePhaseTriplets=[fieldRelativeSpace(:) fieldRelativeTime(:) absSpaceTimePhaseTriplets(:,3)];
 
       currFieldSpaceTimePhaseTriplets=[absSpaceTimePhaseTriplets(:,1) absSpaceTimePhaseTriplets(:,2) absSpaceTimePhaseTriplets(:,3)];


	allCellExcessJoint=allCellExcessJoint+excessJointSpikeDensity;

   numCellsUsed=numCellsUsed+1;

   if(plotIndividualPhasePrecession)
    showInfoPlots=1;
    [IposPhase,Hpos,Hphase,rhoSpear,pSpear,dispXlims,dispYlims]=getDiscreteMutualInfo(linearPos,phasesInRadians,30,30,'position','phase (radians)',showInfoPlots);
      %subplot(2,2,2)
      subplot(3,2,2)
	 plot(currFieldSpaceTimePhaseTriplets(:,1),currFieldSpaceTimePhaseTriplets(:,3)*pi/180,'k.')
     xlim(dispXlims)
     ylim(dispYlims)
     xlabel('position in track')
     ylabel('theta phase (radians')
	title(sprintf('p_{corr}=%0.5f,slope=%.5f deg/field, rho=%.5f',pK,sK*360,rhoK))
       drawnow
       maxFig
	setFigFontTo(18)

	saveas(gcf,sprintf('testIndividualPhasePrecession%d.tif',i))
	close all
	continue
   end
       
   
   if(plotIndividualCells)
       [zs] = threeDPtsToSurface(currFieldSpaceTimePhaseTriplets')
       avgPhasePerSpaceTimeBin=NaN(size(zs));
        for row=1:size(zs,1)
            for col=1:size(zs,2)
                %if(~isempty(zs{row,col}))
                if(length(zs{row,col})>0)
                 avgPhasePerSpaceTimeBin(row,col)=circMeanDeg(zs{row,col});
                end
            end
        end

        plotSpaceTimeRF(spaceBins,timeBins,avgPhasePerSpaceTimeBin,fT)
        saveas(gcf,sprintf('testIndividualFields%d.tif',i))
   end
   
    allFieldSpaceTimePhaseTriplets=[allFieldSpaceTimePhaseTriplets; currFieldSpaceTimePhaseTriplets];
    size(allFieldSpaceTimePhaseTriplets)
end
%%
%figure;
%histogram(estSpeed)

if(combineDirections==1 && runDir==1)
    continue
end

%bins using unique function - round to desired bin
allFieldSpaceTimePhaseTriplets(:,1)=round(allFieldSpaceTimePhaseTriplets(:,1)*numSpaceBins);
allFieldSpaceTimePhaseTriplets(:,2)=round(allFieldSpaceTimePhaseTriplets(:,2)*numTimeBins);

allFieldSpaceTimePhaseTriplets(any(isnan(allFieldSpaceTimePhaseTriplets),2),:)=[];

%figure; omarPcolor((1:numDensityBins)/numDensityBins,(1:numDensityBins)/numDensityBins*2*pi,allCellExcessJoint/numCellsUsed)
[X,Y]=meshgrid((1:numDensityBins)/numDensityBins,(1:numDensityBins)/numDensityBins*2*pi);

figure;surf(X,Y,allCellExcessJoint/numCellsUsed)
xlabel('Position in field')
ylabel('spike phase')
zlabel('Excess joint density')
colormap(jet)
colorbar
title(sprintf('N=%d precessing place cells',numCellsUsed))

save('allFieldSpaceTimePhaseTripletsAndExcessJoint.mat','allFieldSpaceTimePhaseTriplets','allCellExcessJoint','numCellsUsed')
fds
%{
if(trackDir==1)
    %plotPolarSpaceTimeField((1:(numSpaceBins+1))-0.5,(1:(numTimeBins+1))-0.5, allFieldSpaceTimePhaseTriplets);
    plotSpacetimePhaseSurface(spaceBins,timeBins,allFieldSpaceTimePhaseTriplets)
else
    maxDist=max(allFieldSpaceTimePhaseTriplets(:,1));
    %allFieldSpaceTimePhaseTriplets(:,1)=maxDist-allFieldSpaceTimePhaseTriplets(:,1);
    maxTime=max(allFieldSpaceTimePhaseTriplets(:,2));
    %allFieldSpaceTimePhaseTriplets(:,2)=maxTime-allFieldSpaceTimePhaseTriplets(:,2);
    %plotPolarSpaceTimeField(((numSpaceBins+1):-1:1)-0.5,(1:(numTimeBins+1))-0.5, allFieldSpaceTimePhaseTriplets);
%}  
fH=figure;

if(dummyRun)
    minValData=148.3942;
    maxValData=227.1075;
    
    dummyPhasePerPos=linspace(maxValData,minValData,numSpaceBins+1);
    
    smoothPhaseSurface=ones(numSpaceBins+1,numTimeBins+1);
    for si=1:numSpaceBins+1
        smoothPhaseSurface(si,:)=smoothPhaseSurface(si,:)*dummyPhasePerPos(si);
    end
    
    subplot(2,2,1)
        
       spaceBins=1:(numSpaceBins+1);
       timeBins=1:(numTimeBins+1);
       spaceBins=scaledata(spaceBins,0,1);
        timeBins=scaledata(timeBins,0,1);
    
        %surf(X,Y,smoothPhaseSurface')
         omarPcolor(spaceBins,timeBins,smoothPhaseSurface',fH)
         %omarPcolor(timeBins,spaceBins,smoothPhaseSurface,fH)
        colormap(jet)
        colorbar
        
        
        %[X,Y]=meshgrid(spaceBins(1:maxSpaceID)-0.5,timeBins(minTimeID:(end))-0.5);
        [X,Y]=meshgrid(spaceBins,timeBins);
        
        caxis([minPhaseDisp maxPhaseDisp])
       
        %zlim([50 300])
        

            xlabel('Norm distance from left edge of field')
         ylabel('Norm time into field')
        
        %xlim([160 160+310])
        %subplot(2,1,2)
         subplot(2,2,2)
        [~,cpH]=contour(X,Y,smoothPhaseSurface',30)

           
           caxis([minPhaseDisp maxPhaseDisp])
         title('Isophase lines')
         %xlabel('Distance from start wall (relative to field)')
         xlabel('Norm distance from left edge of field')
         ylabel('Norm time into field')
         cpH.LineWidth=5;
        
        colormap(jet)
        colorbar
else
    [smoothPhaseSurface]=plotSpacetimePhaseSurface(spaceBins,timeBins,allFieldSpaceTimePhaseTriplets,fH);
end

%%
subplot(2,2,3)

plotBezierCurvesDiffConcavities
title('Trajectories through spacetime field at different speeds')
colorbar
xlabel('Norm distance from left edge of field')
ylabel('Norm time into field')
         

%%
%figure
subplot(2,2,4)
numConvexities=20;
numConvexities=30;
%numConvexities=50;
convexities=linspace(-1,1,numConvexities);
convexities=linspace(-1,1,numConvexities);
%convexities=linspace(-0.5,0.5,numConvexities);
colormp=copper(numConvexities);

for ci=1:numConvexities
    
    currConvexityFrac=convexities(ci);
%currConvexityFrac=1;
     [curveX, curveY] = getBezierGivenConvexity(-currConvexityFrac);
     curveYScaled=scaledata(curveX,1,size(smoothPhaseSurface,2));
     curveXScaled=scaledata(curveY,1,size(smoothPhaseSurface,1));

    c=improfile(smoothPhaseSurface,curveYScaled,curveXScaled);
    %newXQ=linspace(1,length(curveXScaled),length(c));
    %xPoints=interp1(1:length(curveXScaled),curveXScaled,newXQ);
    
    c=c(~isnan(c));
    newN=100;
    newQ=linspace(1,length(c),newN);
    cResamp=interp1(1:length(c),c,newQ);
    
    hold on
    plot(1:newN,cResamp','Color',colormp(ci,:),'LineWidth',1)
    %plot(xPoints,c,'Color',colormp(ci,:),'LineWidth',3)

end
colormap(gca,copper)
    cb=colorbar
    ylabel(cb,'Convexity of trajectory in spacetime')
    caxis([-1 1])
xlabel('Percent of trajectory traversed')
ylabel('Theta phase')
title('Phase across field trajectory for different speeds')

ylim([minPhaseDisp+30 maxPhaseDisp-20])
%%
%title(sprintf('Average distance-time-phase relationship, direction %d',trackDir))
numDataPts=size(allFieldSpaceTimePhaseTriplets,1);
if(dummyRun==1)
    uberTitle('Control spacetime phase surface')
else
    if(combineDirections)
        uberTitle(sprintf('Spacetime Phase Surface (n=%d place cell-field pairs), both directions',numCellsUsed))
    else
        uberTitle(sprintf('Spacetime Phase Surface (n=%d place cells), direction %d',numCellsUsed,trackDir))
    end
end

%maxFig
setFigFontTo(18)
maxFig

if(dummyRun==1)
    saveas(gcf,sprintf('groupPlaceCellDirSpaceTimePhaseSurfaceControl.tif'))
else  
    saveas(gcf,sprintf('groupPlaceCellSpaceTimePhaseSurface.tif'))
end
%saveas(gcf,sprintf('groupPlaceCellDir%d_SpaceTimePhaseSurface.tif',trackDir))
%{
plotPolarSpaceTimeField(spaceBins,timeBins,allFieldSpaceTimePhaseTriplets)
title(sprintf('Average distance-time-phase relationship (n=%d place cells), direction %d',numCellsUsed,trackDir))
%maxFig
setFigFontTo(18)
saveas(gcf,sprintf('groupPlaceCellDir%d_SpaceTimeAnalysisPolarPlot.tif',trackDir))
%}

end %end direction loop

end %end dummy run loop

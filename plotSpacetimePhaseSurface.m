function [smoothPhaseSurface] = plotSpacetimePhaseSurface(spaceBins,timeBins,allFieldSpaceTimePhaseTriplets,fH)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
setTightSubplots_SpaceTime
    surfMode=0;
    polarMode=1;
    plotTimeAsColor=false;
    %plotTimeAsColor=true;
    %numCellsUsed=;
    %minSampleCount=10;
    minSampleCount=15;
    minSampleCount=1;
    minSampleCount=10;
     minSampleCount=5;
     minSampleCount=1;
    
    minPhaseDisp=40;
    maxPhaseDisp=200;
    
    minPhaseDisp=0;
    maxPhaseDisp=360;
    
%       minPhaseDisp=90;
%    maxPhaseDisp=270;
    
 %    minPhaseDisp=120;
 %   maxPhaseDisp=250;
    
    
    minTime=2;
    maxTime=10;
      %maxTime=8;
      maxTime=15;
    if(~exist('fH','var'))
        fH=figure;
    end
    
    startSpace=spaceBins(1);
    endSpace=spaceBins(2);
    
    startTime=timeBins(1);
    endTime=timeBins(2);
    
        %roundPhases=round(avgPhasePerSpaceTimeBin);
        numPhaseBins=360;
        %numPhaseBins=120;
        numPhaseBins=20;
        numPhaseBins=24;
         %numPhaseBins=10;
        %numPhaseBins=30;
         numPhaseBins=40;
          %numPhaseBins=60;
        %numPhaseBins=15;
        %numPhaseBins=8;
        phaseEdges=linspace(0,360,numPhaseBins+1);
        numSpaceBins=50+1;
        numTimeBins=15+1;
        spaceBins=1:numSpaceBins;
        timeBins=1:numTimeBins;
        numDataPts=length(allFieldSpaceTimePhaseTriplets);
        
        spaceTimeToPhase=NaN(numSpaceBins,numTimeBins,1000);
        spaceTimeToPhaseCount=zeros(numSpaceBins,numTimeBins);
       for i=1:numDataPts
            %[row,col]=ind2sub(size(avgPhasePerSpaceTimeBin),spaceTimesPerPhase(i));
            spaceBin=allFieldSpaceTimePhaseTriplets(i,1)+1;
            timeBin=allFieldSpaceTimePhaseTriplets(i,2)+1;
            currPhase=allFieldSpaceTimePhaseTriplets(i,3);
            
            %if(timeBin==1 || spaceBin==numSpaceBins)
            %    continue
            %end
            

            spaceTimeToPhaseCount(spaceBin,timeBin)=spaceTimeToPhaseCount(spaceBin,timeBin)+1;
            spaceTimeToPhase(spaceBin,timeBin,spaceTimeToPhaseCount(spaceBin,timeBin))=currPhase;
       end

       % countEdges=linspace(0,50,100);
       %figure; histogram(spaceTimeToPhaseCount(:),countEdges)
       %figure; imagesc(spaceTimeToPhaseCount)
        spaceTimeToPhaseCount(spaceTimeToPhaseCount==0)=NaN;
         lowCountSpaceTimePts=spaceTimeToPhaseCount<minSampleCount;
         
        %spaceTimeToPhase=spaceTimeToPhase./spaceTimeToPhaseCount;
        spaceTimeToPhase=circMeanDeg(spaceTimeToPhase,3);
        spaceTimeToPhase(lowCountSpaceTimePts)=NaN;
        
        minDisp=prctile(spaceTimeToPhase(:),1);
        maxDisp=prctile(spaceTimeToPhase(:),99);
        
        %caxis([minDisp maxDisp])
        %caxis([0 1])
        %title(sprintf('N=%d cells',numCellsUsed))
        


    phaseBins=edgesToBins(phaseEdges);

        
        maxSpaceID=numSpaceBins-1;
        minTimeID=2;
        spaceBins=scaledata(spaceBins,0,1);
        timeBins=scaledata(timeBins,0,1);
        
        %[X,Y]=meshgrid(spaceBins(1:maxSpaceID)-0.5,timeBins(minTimeID:(end))-0.5);
        [X,Y]=meshgrid(spaceBins,timeBins);
        smoothSizeSpace=13;
          %smoothSizeSpace=11;
          %smoothSizeSpace=15;
          %smoothSizeSpace=19;
        smoothSizeTime=3;
            smoothSizeTime=5;
        K=1/(smoothSizeSpace*smoothSizeTime)*(ones(smoothSizeSpace,smoothSizeTime));
        
        %spaceTimeToPhase=spaceTimeToPhase(1:maxSpaceID,minTimeID:(end));
        
        %smoothPhaseSurface=conv2(spaceTimeToPhase,K,'same');
        [smoothPhaseSurface] = smoothSpaceTimePhase(spaceTimeToPhase,smoothSizeSpace,smoothSizeTime);
       numTimesResmooth=2;
        numTimesResmooth=1;
          %numTimesResmooth=0;
        for s=1:numTimesResmooth
            %smoothPhaseSurface=conv2(smoothPhaseSurface,K,'same');
            [smoothPhaseSurface] = smoothSpaceTimePhase(smoothPhaseSurface,smoothSizeSpace,smoothSizeTime);
        end

        %surf(X,Y,repeatedPosition)
        %subplot(2,1,1)
         subplot(2,2,1)
        
        
        %surf(X,Y,smoothPhaseSurface')
         omarPcolor(spaceBins,timeBins,smoothPhaseSurface',fH)
         %omarPcolor(timeBins,spaceBins,smoothPhaseSurface,fH)
        colormap(jet)
        colorbar
        
        %minDispSpace=ceil(smoothSizeSpace/2)+numTimesResmooth-1;
        minDispSpace=ceil(smoothSizeSpace/2)+numTimesResmooth;
        minDispTime=ceil(smoothSizeTime/2)+numTimesResmooth+1;
        
        maxDispSpace=maxSpaceID-minDispSpace-numTimesResmooth-2;
         %maxDispSpace=maxSpaceID-minDispSpace-numTimesResmooth;
        maxDispTime=max(timeBins)-minDispTime+1;
         colormap(jet)
         colorbar
        %caxis([0 360])
        
        caxis([minPhaseDisp maxPhaseDisp])
       
        %zlim([50 300])
        
         edgeCutOff=0;
         xlim([edgeCutOff 1-edgeCutOff])
        ylim([edgeCutOff 1-edgeCutOff])
        %xlim([minDispSpace maxDispSpace])
         %  ylim([minDispTime maxDispTime])
            xlabel('Norm distance from left edge of field')
         ylabel('Norm time into field')
        
        %xlim([160 160+310])
        %subplot(2,1,2)
         subplot(2,2,2)
        [~,cpH]=contour(X,Y,smoothPhaseSurface',30)
         %xlim([minDispSpace maxDispSpace])
         %  ylim([minDispTime maxDispTime])
         
        
         xlim([edgeCutOff 1-edgeCutOff])
        ylim([edgeCutOff 1-edgeCutOff])
           
           caxis([minPhaseDisp maxPhaseDisp])
         title('Isophase lines')
         %xlabel('Distance from start wall (relative to field)')
         xlabel('Norm distance from left edge of field')
         ylabel('Norm time into field')
         cpH.LineWidth=5;
        %contour(X,Y,smoothedRepeatedPos,15)
        
        colormap(jet)
        colorbar
        %caxis([0 360])
        %figure; imagesc(spaceTimeToPhaseCount)
        %colormap(hsv)
        %colorbar
        %caxis([minDisp maxDisp])

        %xlim([150 150+360])
        %xlim([150 150+310])
         %%
           
        %[temporalPhaseGradient,spatialPhaseGradient]=imgradientxy(smoothPhaseSurface);%y is rows, x is cols
        subplot(2,2,3)
        %bezier monotonic traversal curves: start at (0,0), end at (1,1)
        
        
        %{
        omarPcolor(spaceBins,timeBins,-spatialPhaseGradient',fH)
         title('Negative phase gradient with distance')
        colorbar
        
         %xlabel('Distance from start wall (relative to field)')
         xlabel('Norm distance from left edge of field')
         ylabel('Norm time into field')
        %}
        
        subplot(2,2,4)
        
        %{
         %omarPcolor(spaceBins,timeBins,-temporalPhaseGradient',fH)
          %title('Negative phase gradient with time')
          [X,Y]=meshgrid(timeBins,spaceBins);
          quiver(X,Y,-spatialPhaseGradient,-temporalPhaseGradient,2.5)
          
          xlim([-Inf Inf])
          ylim([-Inf Inf])
         
         %xlabel('Distance from start wall (relative to field)')
         xlabel('Norm distance from left edge of field')
         ylabel('Norm time into field')
         colorbar
        %}
disp('here')


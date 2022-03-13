function []=plot3DjointDistrSlices(x,y,z,xStr,yStr,zStr)
%UNTITLED14 Summary of this function goes here
%   Detailed explanation goes here
setTightSubplots_SpaceTime
         x(x>1)=NaN;
        y(y>1)=NaN;
        z(z>1)=NaN;
        
        anyNaNrow=(isnan(x) | isnan(y) | isnan(z));
        x(anyNaNrow)=[];
        y(anyNaNrow)=[];
        z(anyNaNrow)=[];
        
        univSmoothParam=0.2;
        numBins=50; %ADJUST IN HISTCN
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %get 2D prior joint distribution
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [pXYn pXYEdges pXYMids xyInBins] = histcn([x(:) y(:)]);
    pXYn=pXYn(1:(end-1),1:(end-1));
    pXY=pXYn/sum(pXYn(:));
    defaultSmoothParamXY=0.0111;
    [pXYSmoothed,smoothParamXY]=smoothn(pXY,univSmoothParam); 
    [pXYSmoothed,smoothParamXY]=smoothn(pXYSmoothed,univSmoothParam); 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %get 3D joint distribution
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [pXYZn pXYZEdges pXYZMids xyzInBins] = histcn([x(:) y(:) z(:)]);
    pXYZn=pXYZn(1:(end-1),1:(end-1),1:(end-1));
    pXYZ=pXYZn/sum(pXYZn(:));
    defaultSmoothParamXYZ=0.0527;
    [pXYZSmoothed,smoothParamXYZ]=smoothn(pXYZ,univSmoothParam); 
    [pXYZSmoothed,smoothParamXYZ]=smoothn(pXYZSmoothed,univSmoothParam); 
    [pXYZSmoothed,smoothParamXYZ]=smoothn(pXYZSmoothed,univSmoothParam); 
    
    pXYSmoothed(pXYSmoothed<0)=0;
    pXYZSmoothed(pXYZSmoothed<0)=0;
    
    %figure
    %imagesc(pXYSmoothed)
    
    %colorbar
    
    fH=figure;
    numPlotsPerRow=5;
    maxDispProb=0.1;
     maxDispProb=0.05;
        %zSlices=linspace(0,max(z(:))-0.01,numBins;
        dispPxyz=pXYZSmoothed;
        %threshDispP=prctile(pXYZ(:),90);
        
        %dispPxyz(dispPxyz<threshDispP)=NaN;
        plotCount=1;
        for xi=1:2:size(dispPxyz,1)
            subplot(numPlotsPerRow,numPlotsPerRow,plotCount)
            currPhaseHist=squeeze(dispPxyz(xi,:,:));
            
            %divide each row by prior probability of current time-dist pair
            for yi=1:numBins
             currPhaseHist(yi,:)=currPhaseHist(yi,:)./pXYSmoothed(xi,yi);
            end
            
            currPhaseHistSmoothed=currPhaseHist;
            currPhaseHistSmoothed=nanGaussSmooth(currPhaseHist);
            
            currPhaseHistSmoothed(currPhaseHistSmoothed>maxDispProb*5)=NaN;
            omarPcolor(pXYZMids{2},pXYZMids{3},currPhaseHistSmoothed',fH)
            %s=slice(pXYZMids{1},pXYZMids{2},pXYZMids{3}, dispPxyz, 0, 0, zSlices(zi));
     
                xlabel(yStr)
                ylabel(zStr)
                %zlabel(zStr)
                colormap(gca,jet)
                cb3=colorbar
                ylabel(cb3,'posterior prob.')
                hold on
                    daspect([1 1 1])
                    %title(sprintf('Distance = %.2f',pXYZMids{1}(xi)))
                     title(sprintf('%s = %.2f',xStr,pXYZMids{1}(xi)))
                    %title('pXYZ')
                  

        %set(s, 'EdgeColor', 'none')
        %view([90 0])
        caxis([0 maxDispProb])
            drawnow
            plotCount=plotCount+1;
            %pause(0.2)
            %alpha 0.3
    
        end
          maxFig

      uberTitle(sprintf('%s vs %s posterior joint distributions, per %s (170 fields)', zStr, lower(yStr), lower(xStr)))
            setFigFontTo(18)
            saveas(gcf,sprintf('%sVs%sPosteriorSliceDistrPer%s.tif',zStr, (yStr), (xStr)))


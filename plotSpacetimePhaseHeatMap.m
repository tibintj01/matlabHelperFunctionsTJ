function [outputDataStruct]=plotSpacetimePhaseHeatMap(data3cols,binSize,summaryFlag,previousOutputDataStruct,ratePerTime,totalFieldCount)



setTightSubplots_SpaceTime

normTimeWithinField=1;
%{
xMin=-0.5
xMax=0.5
tMin=-0.5
tMax=0.5
%}

xMin=0
xMax=1
tMin=0
tMax=1


if(~exist('previousOutputDataStruct','var') || isempty(previousOutputDataStruct))

    inputStruct.data3cols=data3cols;
    %inputStruct.xRad=0.01; %field frac
    %inputStruct.yRad=0.01; %field frac
    if(~exist('binSize'))
        inputStruct.xRad=0.015; %field frac
        inputStruct.yRad=0.015; %field frac
    else
            inputStruct.xRad=binSize; %field frac
        inputStruct.yRad=binSize; %field frac
    end

   %inputStruct.xEdges=linspace(0,1,75+1);
   %inputStruct.yEdges=linspace(0,1,(75+1));
   
   numBins=round(1/binSize);
    %  inputStruct.xEdges=linspace(xMin,xMax,75+1);
   %inputStruct.yEdges=linspace(tMin,tMax,(75+1));
   
   inputStruct.xEdges=linspace(xMin,xMax,numBins+1);
   inputStruct.yEdges=linspace(tMin,tMax,(numBins+1));

   inputStruct.desiredStat='circMean';
   [distBinsPhase,timeBinsPhase,meanPhaseSmooth,dataCountPerBinPhase]=makeHeatMapOf3rdVarUniv(inputStruct);
   
   inputStructRate=inputStruct;
   inputStructRate.desiredStat='mean';
   inputStructRate.data3cols=[data3cols(:,1) data3cols(:,2) ratePerTime(:)];
   [distBinsRate,timeBinsRate,meanRateSmooth,dataCountPerBinRate]=makeHeatMapOf3rdVarUniv(inputStructRate);

   outputDataStruct.distBinsPhase=distBinsPhase;
   outputDataStruct.timeBinsPhase=timeBinsPhase;
   outputDataStruct.meanPhaseSmooth=meanPhaseSmooth;
   outputDataStruct.dataCountPerBinPhase=dataCountPerBinPhase;
   
   outputDataStruct.distBinsRate=distBinsRate;
   outputDataStruct.timeBinsRate=timeBinsRate;
   outputDataStruct.meanRateSmooth=meanRateSmooth;
   outputDataStruct.dataCountPerBinRate=dataCountPerBinRate;

else
   distBinsPhase=previousOutputDataStruct.distBinsPhase;
   timeBinsPhase=previousOutputDataStruct.timeBinsPhase;
   meanPhaseSmooth=previousOutputDataStruct.meanPhaseSmooth;
   dataCountPerBinPhase=previousOutputDataStruct.dataCountPerBinPhase;
   
   distBinsRate=previousOutputDataStruct.distBinsRate;
   timeBinsRate=previousOutputDataStruct.timeBinsRate;
   meanRateSmooth=previousOutputDataStruct.meanRateSmooth;
   dataCountPerBinRate=previousOutputDataStruct.dataCountPerBinRate;
   
   %backwards compatability
   
    outputDataStruct.distBinsPhase=distBinsPhase;
   outputDataStruct.timeBinsPhase=timeBinsPhase;
   outputDataStruct.meanPhaseSmooth=meanPhaseSmooth;
   outputDataStruct.dataCountPerBinPhase=dataCountPerBinPhase;
    outputDataStruct.distBinsRate=distBinsRate;
   outputDataStruct.timeBinsRate=timeBinsRate;
   outputDataStruct.meanRateSmooth=meanRateSmooth;
   outputDataStruct.dataCountPerBinRate=dataCountPerBinRate;

end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %PLOTTING
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(summaryFlag)
   %meanPhaseSmooth(dataCountPerBinPhase<30)=NaN;
      %meanPhaseSmooth(dataCountPerBinPhase<10)=NaN;
      %meanPhaseSmooth(dataCountPerBinPhase<20)=NaN;
      %meanPhaseSmooth(dataCountPerBinPhase<15)=NaN;
      %meanRateSmooth(dataCountPerBinRate<15)=NaN;
       meanPhaseSmooth(dataCountPerBinPhase<10)=NaN;
      meanRateSmooth(dataCountPerBinRate<10)=NaN;
         %meanPhaseSmooth(dataCountPerBinPhase<30)=NaN;
      %meanRateSmooth(dataCountPerBinRate<30)=NaN;
    end
  
   
   fH=figure
       if(summaryFlag)
           %{
        meanPhaseSmoothed=nanGaussSmooth(meanPhaseSmooth);
         meanPhaseSmoothed=nanGaussSmooth(meanPhaseSmoothed);
        meanRateSmoothed=nanGaussSmooth(meanRateSmooth);
          meanRateSmoothed=nanGaussSmooth(meanRateSmoothed);
           %}
          
          meanPhaseSmoothed=nanGaussSmoothCentral(meanPhaseSmooth);
         %meanPhaseSmoothed=nanGaussSmooth(meanPhaseSmoothed);
        %meanRateSmoothed=nanGaussSmoothCentral(meanRateSmooth);
          %meanRateSmoothed=nanGaussSmooth(meanRateSmoothed);
          meanRateSmoothed=nanGaussSmoothRate(meanRateSmooth);

       else
           meanPhaseSmoothed=nanGaussSmoothIndividual(meanPhaseSmooth);
         meanPhaseSmoothed=nanGaussSmoothIndividual(meanPhaseSmoothed);
        meanRateSmoothed=nanGaussSmoothIndividual(meanRateSmooth);
          meanRateSmoothed=nanGaussSmoothIndividual(meanRateSmoothed);
       end
         
       if(summaryFlag)
        subplot(2,2,1)
       else
          subplot(1,2,1)
       end
%omarPcolor(distBinsPhase,timeBinsPhase,meanPhaseSmooth',fH)
omarPcolor(distBinsPhase,timeBinsPhase,meanPhaseSmoothed',fH)
        hold on
     xlabel('Distance in field (fraction)')
         if(normTimeWithinField)
             ylabel('Time in field (fraction)')
         else
             ylabel('Time (sec)')
         end
        colormap(gca,jet)
        cb=colorbar
        ylabel(cb,'Circ mean theta phase (deg)')
   box off

   %caxis([0 360])
   caxis([60 330])
    title({'Theta phase vs space-time'})
    
      xlim([0 1])
         ylim([0 1])
     %plot(xlim,[1 1],'k--','LineWidth',3)
        % plot([1 1],ylim, 'k--','LineWidth',3)

       
            daspect([1 1 1])
            %axis square
   
        if(summaryFlag)
        subplot(2,2,2)
       else
          subplot(1,2,2)
       end
omarPcolor(distBinsRate,timeBinsRate,meanRateSmoothed',fH)
        hold on
     xlabel('Distance in field (fraction)')
         if(normTimeWithinField)
             ylabel('Time in field (fraction)')
         else
             ylabel('Time (sec)')
         end
        colormap(gca,parula)
        cb5=colorbar
        ylabel(cb5,'mean spike rate (Hz)')
   box off
   daspect([1 1 1])
    title({'Rate vs space-time'})
    xlim([0 1])
            ylim([0 1])
            maxRateDisp=prctile(meanRateSmoothed(:),99);
            try
             caxis([0 maxRateDisp])
            end
            
            
 %plot(xlim,[1 1],'k--','LineWidth',3)
  %       plot([1 1],ylim, 'k--','LineWidth',3)
         
   
         
         if(summaryFlag)
             if(exist('totalFieldCount','var'))
                uberTitle(sprintf('Space-time field population mean (n=%d fields)',totalFieldCount))
             else
                   uberTitle({'Space-time field population mean'})
             end
             subplot(2,2,3)
             [X,T]=meshgrid(distBinsPhase,timeBinsPhase);
             %[C,hD]=contour(X,T,meanPhaseSmoothed',30);
               [C,hD]=contour(X,T,meanPhaseSmoothed',40);

           hD.LineWidth=2;
           colormap(gca,jet)
           cb=colorbar
           hold on
            daspect([1 1 1])
         xlim([0 1])
                ylim([0 1])

                title('iso-phase lines')
                 xlabel('Distance in field (frac)')
               ylabel('Time in field (frac)')
               ylabel(cb,'Circ mean theta phase (deg)')
               %caxis([0 360])
                caxis([60 330])


               subplot(2,2,4)
               [Xr,Tr]=meshgrid(distBinsRate,timeBinsRate);
                     [C,hD]=contour(Xr,Tr,meanRateSmoothed',30);

               hD.LineWidth=2;
               colormap(gca,parula)
               cb=colorbar
               hold on
               daspect([1 1 1])
             xlim([0 1])
                        ylim([0 1])

                         xlabel('Distance in field (frac)')
               ylabel('Time in field (frac)')
               ylabel(cb,'mean spike rate (Hz)')

                title('iso-rate lines')
                maxFig
                setFigFontTo(24)
                 saveas(gcf,'SpaceTimePhaseRateWithIsoLines.tif')
         end
              
 %setFigFontTo(18)
    %maxFig
    
   if(summaryFlag)
       
       commonCaxis=caxis
         %subplot(1,2,2)
	figure
	subplot(2,2,2)
         meanPhaseSmoothed=nanGaussSmooth(meanPhaseSmooth);
         meanPhaseSmoothed=nanGaussSmooth(meanPhaseSmoothed);
         %meanPhaseSmoothed=nanGaussSmooth(meanPhaseSmoothed);
   [X,T]=meshgrid(distBinsPhase,timeBinsPhase);
   %[C,hD]=contour(X,T,meanPhaseSmoothed'/360,30);
    [C,hD]=contour(X,T,meanPhaseSmoothed',30);
   
   hD.LineWidth=2;
   colormap(gca,jet)
   cb=colorbar
   hold on
   %plot([-0.5 0.5],[-0.5 0.5],'k--','LineWidth',3)
   
   
   adjX=0.02;
     %plot([-0.12 -0.12]+adjX,ylim,'k--','LineWidth',3)
      % plot(xlim, [-0.22 -0.22],'k--','LineWidth',3)
       
       ptOnLine1=[-0.12+adjX; -0.22];
       ptOnLine2=[.32; .173333333];
       
       diffVector=ptOnLine2-ptOnLine1;
       ptOnLine2=ptOnLine2+diffVector*1.1;
       lineSlope=diffVector(2)/diffVector(1);
       
       %plot([ptOnLine1(1) ptOnLine2(1)], [ptOnLine1(2) ptOnLine2(2)],'k--','LineWidth',3)
       
       adjXLine=0.05;
       plot([xMin xMax],[tMin tMax],'k--','LineWidth',3)
       %plot([-0.5 0.5],[-0.5 0.35],'k--','LineWidth',3)
       %plot([-0.5+adjXLine 0.5+adjXLine],[-0.5 0.5],'--','Color',getGrayRGB,'LineWidth',3)
      
   %caxis([0 1])
   xlabel('Distance in field (frac)')
   ylabel('Time in field (frac)')
   %ylabel(cb,'Avg phase in spacetime (cycle frac)')
   ylabel(cb,'Avg phase in spacetime (deg)')
   title('isophase lines')
   %title(sprintf('%d place cell grouped spacetime response',numFiles))
   %xlim([0 1])
    %ylim([0 1])
    %caxis(commonCaxis)
   hold on
%plot([0 0.7],[0 1],'k--','LineWidth',3)
%plot([0 1],[0 0.7],'k--','LineWidth',3)
%plot([0 1],[0 1],'k--','LineWidth',2)
box off
   %daspect([1 1 1])
   axis square
    

subplot(2,2,1)
   surf(X,T,meanPhaseSmoothed')

   colormap(gca,jet)
   cb=colorbar
   %caxis([0 1])
   %xlabel('Distance in field (frac)')
   ylabel('Time in field (frac)')
   ylabel(cb,'Avg phase in spacetime (deg)')
   %title('isophase lines')
   %title(sprintf('%d place cell grouped spacetime response',numFiles))
   %xlim([0 1])
    %ylim([0 1])
   hold on
%plot([0 0.7],[0 1],'k--','LineWidth',3)
%plot([0 1],[0 0.7],'k--','LineWidth',3)
%plot([0 1],[0 1],'k--','LineWidth',2)
box off

stretchFact=1.15*sqrt(2);
xPts=X(:);
tPts=T(:);
zPts=zeros(size(tPts));
originalSize=size(X);

spacetimePts=[xPts*stretchFact tPts zPts];
spacetimePtsRotated=rotz(45)*spacetimePts';

%rotatedBasis=[sqrt(2)/2 -sqrt(2)/2 0; sqrt(2)/2 sqrt(2)/2 0; 0 0 1];
%backToOurBasis=inv(rotatedBasis);


%spacetimePtsStretched=backToOurBasis*[1 0 0; 0 stretchFact 0; 0 0 1]*rotatedBasis*spacetimePts';
%// Then after transformation
Xnew = reshape(spacetimePtsRotated(1,:), originalSize);
Tnew = reshape(spacetimePtsRotated(2,:), originalSize);
%Znew = reshape(spacetimePtsRotated(:,3), original_size);
daspect([1 1 360])
 view(0,90)
 %{
subplot(2,2,3)
maxVal=max(max(Tnew.^2-Xnew.^2));

conePlotVals=X.^2-T.^2;
conePlotVals(conePlotVals<0)=NaN;
domainIdxes=X/stretchFact>=xMin & X/stretchFact<=xMax & T>=tMin & T<=tMax;
conePlotVals(~domainIdxes)=NaN;
conePlotVals=scaledata(conePlotVals,0,360)

   %surf(X,T,360-(Tnew.^2-Xnew.^2)/maxVal*360)
   surf(Xnew/stretchFact*sqrt(2),Tnew/stretchFact*sqrt(2),360-conePlotVals)
hold on
   surf(Tnew/stretchFact*sqrt(2),Xnew/stretchFact*sqrt(2),360-conePlotVals)

   colormap(gca,jet)
   cb=colorbar
   caxis([120 360])
   %xlabel('Distance in field (frac)')
   ylabel('Time in field (frac)')
   ylabel(cb,'X^2-T^2')
   %title('isophase lines')
   %title(sprintf('%d place cell grouped spacetime response',numFiles))
   xlim([xMin xMax])
    ylim([tMin tMax])
   hold on
   title('Transformed x^2-t^2 cone')
   
%plot([0 0.7],[0 1],'k--','LineWidth',3)
%plot([0 1],[0 0.7],'k--','LineWidth',3)
%plot([0 1],[0 1],'k--','LineWidth',2)
box off
daspect([1 1 360])
   %daspect([1 1 1])
   
    %axis([0 1000 -0.001 110 0 110])    
 view(0,90)
 %}
   setFigFontTo(28)
    maxFig
    saveas(gcf,'SpaceTimePhasePrecessionAvg.tif')
   
   else
	%{
       subplot(1,2,2); 
       plot(data3cols(:,1),data3cols(:,3),'k.')
       ylim([0 360])
       xlim([xMin xMax])
	%}
        maxFig
        subplot(1,2,1); 
       
   end
   
   setFigFontTo(32)

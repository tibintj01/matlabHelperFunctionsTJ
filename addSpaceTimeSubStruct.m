close all; clear all; clc
dataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles';
%xxxxxsetTightSubplots_SpaceTime

filePaths=getFilePathsRegex(dataDir,'*mat');

saveFileName='spaceTimeDataInTheta.mat';
showPhaseVsTime=0;
plotAll3DpointsScatter=0;

if(~exist(saveFileName,'file'))

startFile=1;
useRefThetaCh=1;
justFirstSpikes=0;

minTimeGain=-50;
minSpaceGain=-50;

allX=[];
allTimeTheta=[];
allPhase=[];

plotCount=0;

badPlotIdxes=[56 60 62 64 65 66 77 81 84 87 88 96 104 113]; 
allTimePhaseSlopes=[];
allSpacePhaseSlopes=[];

for fi=startFile:length(filePaths)
    fi
    currFilePath=filePaths{fi};
    currFileName=getFileNameFromPath(currFilePath);
    fileBaseName=currFileName(1:(end-4));
    
    [currSessionNum,currSessName]=getSessionNum(fileBaseName);
    [leftwardRefUnitData,rightwardRefUnitData]=getRefUnitData(currSessName);
    
    data=load(currFilePath)

    
    for di=1:2
        if(di==2)
            currFieldDirection='leftward';

        elseif(di==1)
            if(strContainsCircSessionName(currFileName))
                continue %circular mazes always leftward laps
            end
            currFieldDirection='rightward';
        end
        
        if(~isfield(data.thetaBasedTimeVsPhaseInfo,currFieldDirection))
            continue
        end
        
        numFields=length(fieldnames(data.thetaBasedTimeVsPhaseInfo.(currFieldDirection)));
        
        heatMapInput.xRad=0.05; %m
        %heatMapInput.yRad=0.1; %theta cycles
         heatMapInput.yRad=0.05; %theta cycles
        numBins=30;
        heatMapInput.xEdges=linspace(0,1,numBins+1);
        heatMapInput.yEdges=linspace(0,1,numBins+1);
        heatMapInput.desiredStat='circMean';
        
        for fii=1:numFields
            currFieldSpaceTimeThetaData=data.thetaBasedTimeVsPhaseInfo.(currFieldDirection).(sprintf('field%d',fii));
            xInField=currFieldSpaceTimeThetaData.allLapInFieldSpikeDistsFieldFrac;
            durationsPerLap=currFieldSpaceTimeThetaData.allLapInFieldNumCycles;
            timeToCycleFact=1/nanmean(currFieldSpaceTimeThetaData.allLapInFieldMeanCycleDurs);
            if(isempty(xInField))  
                continue
            end
            %xInField=scaledata(xInField,0,1);
            xInField=xInField-min(xInField);
            xInField=xInField/max(xInField);
            %xInField=normRangeByPercentile(xInField,10,90);
            
            
            yInField=(currFieldSpaceTimeThetaData.allLapSpikeCycleIDsFromFieldEntry);
            %maxTime=prctile(yInField,80);
            
            %maxTime=prctile(yInField,90);
            %maxTime=prctile(durationsPerLap,75);
            %normDuration=nanmedian(durationsPerLap);
            
            manualTimeFieldEndName=sprintf('manualTimeField%dEnd%sSec',fii,currFieldDirection);
            if(~isfield(data,manualTimeFieldEndName))
                continue
            end
            manualMaxFieldTime=data.(manualTimeFieldEndName);
            
            
       
            
            normDurationInCycles=manualMaxFieldTime*timeToCycleFact;
            
            figure(121)
            subplot(14,14,plotCount+1)
            histogram(durationsPerLap,30)
            hold on
            plot([normDurationInCycles normDurationInCycles],ylim,'k--','LineWidth',3)
            
            %maxTime=max(yInField(:));
            yInField=yInField/normDurationInCycles;
               %yInField=yInField/15;
            %yInField=normRangeByPercentile(yInField,10,90);
            
            zInField=currFieldSpaceTimeThetaData.allLapInFieldSpikePhases;
            
            nonNanIdxes=~isnan(xInField) & ~isnan(yInField) & ~isnan(zInField);
            %distTimeAvg=(xInField+yInField)/2;
            %distTimeAvg=1-(xInField+yInField)/2;
            conePredictionPerPt=NaN(size(xInField));
            spaceOnlyPredictionPerPt=NaN(size(xInField));
            
            for xii=1:length(xInField)
                conePredictionPerPt(xii)=getConeDegPredForXandTvals(xInField(xii),yInField(xii));
                spaceOnlyPredictionPerPt(xii)=getSpaceOnlyDegPredForXandTvals(xInField(xii),yInField(xii));
            end
            %distTimeAvg=(xInField^2-yInField^2);
            
            rMatDistPhase=corrcoef(xInField(nonNanIdxes),zInField(nonNanIdxes));
            rMatTimePhase=corrcoef(yInField(nonNanIdxes),zInField(nonNanIdxes));
            %rMatDistTimePhase=corrcoef(conePredictionPerPt(nonNanIdxes),zInField(nonNanIdxes));
            
            rDistPhase=rMatDistPhase(1,2);
            rTimePhase=rMatTimePhase(1,2);
            %rDistTimePhase=rMatDistTimePhase(1,2);
            
            circRDistPhase=kempter_lincirc(xInField(nonNanIdxes),ang2rad(zInField(nonNanIdxes)));
            circRTimePhase=kempter_lincirc(yInField(nonNanIdxes),ang2rad(zInField(nonNanIdxes)));
            %circRDistTimePhase=kempter_lincirc(conePredictionPerPt(nonNanIdxes),ang2rad(zInField(nonNanIdxes)));

            
            
            fH=figure(122)
            subplot(2,2,1)
            plot(xInField,zInField,'k.')
             xlim([0 1])
            ylim([0 360])
            xlabel('Dist in field (frac)')
            ylabel('Reference theta phase (deg)')
             daspect([1 360 1])
             
             %title(sprintf('R=%.2f',rDistPhase))
             title(sprintf('circR=%.2f',circRDistPhase))
             
             
            
            subplot(2,2,4)
            if(showPhaseVsTime)
                plot(yInField,zInField,'k.')
                xlim([0 1])
                ylim([0 360])
                xlabel('Time in field (frac)')
                ylabel('Reference theta phase (deg)')
                daspect([1 360 1])

                 %title(sprintf('R=%.2f',rTimePhase))
                 title(sprintf('circR=%.2f',circRTimePhase))
            else
                
                
                heatMapInput.data3cols=[xInField(:), yInField(:), spaceOnlyPredictionPerPt(:)];
                [xCenters,yCenters,spaceTimePhaseSpaceOnlyPredicted,dataCountPerBin] =makeHeatMapOf3rdVarUniv(heatMapInput);

                omarPcolor(yCenters,xCenters,spaceTimePhaseSpaceOnlyPredicted',fH)
                colormap(jet)
                cb=colorbar
                caxis([0 360])

                daspect([1 1 1])
                xlabel('Dist in field (frac)')
                ylabel('Time in field (frac)')
                ylabel(cb,'Space only prediction(deg)')
          
                title('Space only model phase prediction')
                setFigFontTo(24)
            end
             %{
            subplot(2,2,3)
            plot(conePredictionPerPt,zInField,'k.')
            xlim([0 1])
            ylim([0 360])
            xlabel('spacetime interval')
            ylabel('Reference theta phase (deg)')
            daspect([1 360 1])
             %title(sprintf('R=%.2f',rDistTimePhase))
             title(sprintf('circR=%.2f',circRDistTimePhase))
             %}
             
            %{
            scatter(xInField(:),yInField(:),20,zInField(:))
            hold on
            colormap(jet)
            cb2=colorbar
              caxis([0 360])
            xlim([0 1])
            ylim([0 1])
            daspect([1 1 1])
             xlabel('Dist in field (frac)')
           ylabel('Time in field (frac)')
            ylabel(cb2,'Reference theta phase (deg)')
            %}
            
            subplot(2,2,3)
            
            heatMapInput.data3cols=[xInField(:), yInField(:), conePredictionPerPt(:)];
            [xCenters,yCenters,spaceTimePhaseConePredicted,dataCountPerBin] =makeHeatMapOf3rdVarUniv(heatMapInput);
            
            omarPcolor(yCenters,xCenters,spaceTimePhaseConePredicted',fH)
            colormap(jet)
            cb=colorbar
            caxis([0 360])
            
            daspect([1 1 1])
            xlabel('Dist in field (frac)')
            ylabel('Time in field (frac)')
            ylabel(cb,'Spacetime Cone prediction(deg)')
            title('Spacetime cone model phase prediction')
            setFigFontTo(24)
        
            subplot(2,2,2)
            heatMapInput.data3cols=[xInField(:), yInField(:), zInField(:)];

            [xCenters,yCenters,spaceTimePhaseMeansSmooth,dataCountPerBin] =makeHeatMapOf3rdVarUniv(heatMapInput);
            [X,Y]=meshgrid(xCenters,yCenters);
            nonNaNIdxes=~isnan(X(:)) & ~isnan(Y(:)) & ~isnan(spaceTimePhaseMeansSmooth(:));
            
            points3DspaceNoNan=[X(nonNaNIdxes), Y(nonNaNIdxes), spaceTimePhaseMeansSmooth(nonNaNIdxes)];
            
            %[bestFitPlaneNormalVector,V,pointOnBestFitPlane] = affine_fit(points3DspaceNoNan);
            predictorMatrix=[ones(size(spaceTimePhaseMeansSmooth(nonNaNIdxes))),points3DspaceNoNan(:,1),points3DspaceNoNan(:,2)];
            bestFitLineCoeff=regress(points3DspaceNoNan(:,3),predictorMatrix);
            numPtsUsedToEstSlope=length(points3DspaceNoNan(:,1));
            
            if(bestFitLineCoeff(2)>minSpaceGain || bestFitLineCoeff(3) > minTimeGain)
                continue
            end

            
            spaceTimePhaseMeansSmooth=nanGaussSmooth(spaceTimePhaseMeansSmooth);
            
            omarPcolor(yCenters,xCenters,spaceTimePhaseMeansSmooth',fH)
            colormap(jet)
            cb=colorbar
            caxis([0 360])
            
            daspect([1 1 1])
            xlabel('Dist in field (frac)')
            ylabel('Time in field (frac)')
            ylabel(cb,'Reference theta phase (deg)')
            title('Spacetime circ avg. phase')
    
            
            plotCount=plotCount+1;
            uberTitle(removeUnderscores(sprintf('%s, num %d',fileBaseName,plotCount)))
            
            subplot(2,2,1)
            
            %plotPlaneSurface(bestFitPlaneNormalVector,pointOnBestFitPlane)
            scatter3(X(nonNaNIdxes), Y(nonNaNIdxes), spaceTimePhaseMeansSmooth(nonNaNIdxes))
            
            hold on
            startPoint=[0,0,bestFitLineCoeff(1)];
            endPoint=[1,1,bestFitLineCoeff(1)+bestFitLineCoeff(2)+bestFitLineCoeff(3)];
            plot3([startPoint(1) endPoint(1)],[startPoint(2) endPoint(2)],[startPoint(3) endPoint(3)],'k-','LineWidth',5)
                    setFigFontTo(14)
                    zlim([0 360])
                    xlim([0 1])
                    ylim([0 1])
                    xlabel('Dist in field (frac)')
            ylabel('Time in field (frac)')
            zlabel('Theta phase (deg)')
                    
                    daspect([1 1 360])
                    view(-28.3000, 10.8000)
                    currTimePhaseSlope=(bestFitLineCoeff(3));
                    currSpacePhaseSlope=(bestFitLineCoeff(2));
                    title({'Spacetime vs phase, 3D linear fit',sprintf('Time-phase slope: %d deg/field, Space-phase slope: %d deg/field', round(bestFitLineCoeff(3)), round(bestFitLineCoeff(2)))})
            maxFig
            
            
            singleFieldConeMatchStats.(sprintf('spaceTimePhaseMeansSmoothDir%dField%d',di,fii))=spaceTimePhaseMeansSmooth;
            singleFieldConeMatchStats.(sprintf('spaceTimePhaseConePredictedDir%dField%d',di,fii))=spaceTimePhaseConePredicted;
             singleFieldConeMatchStats.(sprintf('spaceTimePhaseSpaceOnlyPredictedDir%dField%d',di,fii))=spaceTimePhaseSpaceOnlyPredicted;
             singleFieldConeMatchStats.(sprintf('bestFitLineCoeffDir%dField%d',di,fii))=bestFitLineCoeff;
              %singleFieldConeMatchStats.(sprintf('bestFitPlaneNormalVectorDir%dField%d',di,fii))=bestFitPlaneNormalVector;
             %singleFieldConeMatchStats.(sprintf('pointOnBestFitPlaneDir%dField%d',di,fii))=pointOnBestFitPlane;
             
            singleFieldConeMatchStats.xCenters=xCenters;
            singleFieldConeMatchStats.yCenters=yCenters;
            
            minNumPts3DLineFit=60;
          if(numPtsUsedToEstSlope>minNumPts3DLineFit)
            allTimePhaseSlopes=[allTimePhaseSlopes; currTimePhaseSlope];
            allSpacePhaseSlopes=[allSpacePhaseSlopes; currSpacePhaseSlope];
          end
            
            saveas(gcf,sprintf('spaceTimePhaseInThetaCycles%s_%s_field%d.png',fileBaseName,currFieldDirection,fii))
            
            close(fH)
            if(plotAll3DpointsScatter)
             figure(81)
            
           
            scatter3(points3DspaceNoNan(:,1),points3DspaceNoNan(:,2),points3DspaceNoNan(:,3),3,points3DspaceNoNan(:,3))
            colormap(gca,jet)
             hold on
            xlim([0 1])
            ylim([0 1])
            zlim([0 360])
            end
            
            if(~ismember(plotCount,badPlotIdxes))
                allX=[allX; xInField(:)];
                allTimeTheta=[allTimeTheta; yInField(:)];
                allPhase=[allPhase; zInField(:)];
            end
            
        end
        
        
    end
      save(currFilePath, 'singleFieldConeMatchStats','-append')
end

    save(saveFileName)
else
    load(saveFileName)
end
    
%%
figure; plot(allSpacePhaseSlopes,allTimePhaseSlopes,'k.','MarkerSize',30); 
xlim([-360 120]);  ylim([-360 120]); hold on; plot([-360 120], [-360 120],'k--'); 
daspect([1 1 1])
xlabel('Space-phase 3D linear fit slope (deg/field)')
ylabel('Time-phase 3D linear fit slope (deg/field)')
title('Spacetime vs theta phase; spatial vs temporal gain')
setFigFontTo(18)
saveas(gcf,'spacetimeVsThetaPhaseSpatialVsTemporalGain.png')

fHMaster=figure
heatMapInput.yRad=0.05;
 heatMapInput.data3cols=[allX(:), allTimeTheta(:), allPhase(:)];
 %heatMapInput.data3cols=[allX(:), allTimeTheta(:), 360*(1-(allTimeTheta(:) + allX(:))/2) ];
 
 
            [xCenters,yCenters,spaceTimePhaseMeans,dataCountPerBin] =makeHeatMapOf3rdVarUniv(heatMapInput);
        minNumPts=200;
        %minNumPts=300;
                minNumPts=25;
                %minNumPts=150;
        
        goodIdxes=dataCountPerBin>minNumPts;
        
        spaceTimePhaseMeansSmooth=nanGaussSmooth(spaceTimePhaseMeans);
        
        spaceTimePhaseMeansSmooth(~goodIdxes)=NaN;
            subplot(2,1,1)
            omarPcolor(yCenters,xCenters,spaceTimePhaseMeansSmooth',fHMaster)
            colormap(gca,jet)
            cb2=colorbar
            %caxis([50 360])
              caxis([100 330])
                  xlabel('Dist in field (frac)')
           ylabel('Time in field (frac)')
            ylabel(cb2,'Reference theta phase (deg)')
            
            setFigFontTo(24)
            
            daspect([1 1 1])
            subplot(2,1,2)
            
                  omarPcolor(yCenters,xCenters,dataCountPerBin',fHMaster)
            colormap(gca,parula)
            cb3=colorbar
            %caxis([50 360])
            
            daspect([1 1 1])
            
            xlabel('Dist in field (frac)')
           ylabel('Time in field (frac)')
            ylabel(cb3,'Data point count')
            uberTitle('spacetimeByThetaCycleSummary')
            setFigFontTo(24)
            
            
            
            saveas(gcf,'spacetimeByThetaCycleSummary.tif')

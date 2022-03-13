close all; clear all; clc
dataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles';
%xxxxxsetTightSubplots_SpaceTime

filePaths=getFilePathsRegex(dataDir,'*mat');

saveFileName='spaceTimeDataInTheta.mat';
showPhaseVsTime=0;

startFile=1;
useRefThetaCh=1;
justFirstSpikes=0;

allX=[];
allTimeTheta=[];
allPhase=[];

plotCount=0;

onlyHighSpeedComparison=0;

badPlotIdxes=[56 60 62 64 65 66 77 81 84 87 88 96 104 113]; 

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
            if(~isfield(data.singleFieldConeMatchStats,sprintf('spaceTimePhaseMeansSmoothDir%dField%d',di,fii)))
                continue
            end
            
            xCenters=data.singleFieldConeMatchStats.xCenters;
            yCenters=data.singleFieldConeMatchStats.yCenters;
            
            [X,Y]=meshgrid(xCenters,yCenters)
            
            currFieldSpaceTimePhaseMeans=data.singleFieldConeMatchStats.(sprintf('spaceTimePhaseMeansSmoothDir%dField%d',di,fii));
            
            currFieldSpaceTimePhaseConePredicted=data.singleFieldConeMatchStats.(sprintf('spaceTimePhaseConePredictedDir%dField%d',di,fii));
            currFieldSpaceTimePhaseSpaceOnlyPredicted=data.singleFieldConeMatchStats.(sprintf('spaceTimePhaseSpaceOnlyPredictedDir%dField%d',di,fii));
        
            goodIdxes=~isnan(currFieldSpaceTimePhaseMeans) & ~isnan(currFieldSpaceTimePhaseConePredicted) & ~isnan(currFieldSpaceTimePhaseSpaceOnlyPredicted);
        
            highSpeedIdxes=X<=Y;
            
            if(onlyHighSpeedComparison)
                goodIdxes=goodIdxes & highSpeedIdxes;
            end
          
              if(isempty(currFieldSpaceTimePhaseMeans(goodIdxes)))
                continue
              end
            
            
            dataConeCorrMat=corrcoef(currFieldSpaceTimePhaseMeans(goodIdxes),currFieldSpaceTimePhaseConePredicted(goodIdxes));
            dataSpaceOnlyCorrMat=corrcoef(currFieldSpaceTimePhaseMeans(goodIdxes),currFieldSpaceTimePhaseSpaceOnlyPredicted(goodIdxes));
            coneSpaceOnlyCorrMat=corrcoef(currFieldSpaceTimePhaseConePredicted(goodIdxes),currFieldSpaceTimePhaseSpaceOnlyPredicted(goodIdxes));
            spaceOnlyConeCorrMat=corrcoef(currFieldSpaceTimePhaseSpaceOnlyPredicted(goodIdxes),currFieldSpaceTimePhaseConePredicted(goodIdxes));
            
            
          
            %{
             dataConeCorrMat=corr([currFieldSpaceTimePhaseMeans(goodIdxes),currFieldSpaceTimePhaseConePredicted(goodIdxes)],'Type','Spearman');
            dataSpaceOnlyCorrMat=corr([currFieldSpaceTimePhaseMeans(goodIdxes),currFieldSpaceTimePhaseSpaceOnlyPredicted(goodIdxes)],'Type','Spearman');
            coneSpaceOnlyCorrMat=corr([currFieldSpaceTimePhaseConePredicted(goodIdxes),currFieldSpaceTimePhaseSpaceOnlyPredicted(goodIdxes)],'Type','Spearman');
            spaceOnlyConeCorrMat=corr([currFieldSpaceTimePhaseSpaceOnlyPredicted(goodIdxes),currFieldSpaceTimePhaseConePredicted(goodIdxes)],'Type','Spearman');

            %}
            dataConeCorrR=dataConeCorrMat(2,1);
            dataSpaceOnlyCorrR=dataSpaceOnlyCorrMat(2,1);
            coneSpaceOnlyCorrR=coneSpaceOnlyCorrMat(2,1);
            spaceOnlyConeCorrR=spaceOnlyConeCorrMat(2,1);
            
            currFieldSpaceTimePhaseMeans(~goodIdxes)=NaN;
            currFieldSpaceTimePhaseConePredicted(~goodIdxes)=NaN;
            currFieldSpaceTimePhaseSpaceOnlyPredicted(~goodIdxes)=NaN;
            
            fH=figure(12)
            titleStrs={'Spacetime phase values',sprintf('Cone model phase prediction, R=%.3f',dataConeCorrR),sprintf('Space only model phase prediction, R=%.3f',dataSpaceOnlyCorrR)};
            for pi=1:3
                subplot(2,2,pi)
                if(pi==1)
                    omarPcolor(yCenters,xCenters,currFieldSpaceTimePhaseMeans',fH)
                elseif(pi==2)
                    omarPcolor(yCenters,xCenters,currFieldSpaceTimePhaseConePredicted',fH)
                else
                    omarPcolor(yCenters,xCenters,currFieldSpaceTimePhaseSpaceOnlyPredicted',fH)
                end
                
                colormap(jet)
                cb=colorbar
                caxis([0 360])

                daspect([1 1 1])
                xlabel('Dist in field (frac)')
                ylabel('Time in field (frac)')
                ylabel(cb,'Theta phase (deg)')
          
                title(titleStrs{pi})
                
            end %model loop
            uberTitle({removeUnderscores(sprintf('%s',fileBaseName)),sprintf('cone-space prediction R=%.3f',coneSpaceOnlyCorrR)})
            
            %s=subplot(2,2,4)
           %title({sprintf('cone-space prediction R=%.3f',coneSpaceOnlyCorrR), sprintf('space-cone prediction R=%.3f',spaceOnlyConeCorrR)})
            setFigFontTo(18)
            maxFig
           
            saveas(gcf,sprintf('coneVsSpacePredictions_spaceTimePhaseInThetaCycles%s_%s_field%d.png',fileBaseName,currFieldDirection,fii))

            close all
        end %field loop
        
    end %direction loop
    
end %cell loop
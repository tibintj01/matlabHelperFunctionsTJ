close all; clear all; clc

processedDataDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData';
timeBoundSaveDir='./hc3Images/timeVsPhaseWithBound';

touchDir(timeBoundSaveDir)

unitDataDir=fullfile(processedDataDir,'unitSpikeInfoStructs');

unitFilePaths=getFilePathsRegex(unitDataDir,'*mat');


plotIndFields=1;
plotIndFields=0;

totalFieldCount=0;
%setTightSubplots
maximizeInitializePhase=0;
%maximizeInitializePhase=1;
numTimeBins=120;
numTimeBins=60;
%numTimeBins=40;

minNumSpikes=500;

minNumSpikes=1000;

%minNumSpikes=900;
minNumSpikes=0;

minRunSpeedDisp=0.05;%m/s
maxRunSpeedDisp=0.8; %m/s

plotByLap=0; 
plotBySpeed=0;
plotByDist=1;

startFi=1;

allUnitSpikeTimeFieldFracs=[];
allUnitSpikeTimeFieldPhases=[];

figure

for fi=startFi:length(unitFilePaths)
    
    if(mod(fi,50)==0)
        disp(fi)
    end
    currUnitFilePath=unitFilePaths{fi};
    currUnitStruct=load(currUnitFilePath);
    
    if(isfield(currUnitStruct,'cycleTimeFieldBoundPerDirPerField'))
        currUnitData=currUnitStruct.unitInfo;

        currUnitFileName=getFileNameFromPath(currUnitFilePath);
        
        cycleTimeFieldBoundPerDirPerField=currUnitStruct.cycleTimeFieldBoundPerDirPerField;

        numRightwardFields=length(currUnitStruct.rightwardFieldStartEndM)/2;
        numLeftwardFields=length(currUnitStruct.leftwardFieldStartEndM)/2;

        maxNumFieldsPerDir=max([numRightwardFields numLeftwardFields]);

        %if(redoBounds)
            rightwardTimeBounds=NaN(maxNumFieldsPerDir,1);
        %end
        
        for di=1:2
            if(di==1)
                currFieldDirStr='rightward';
                numFields=numRightwardFields;
            else
                currFieldDirStr='leftward';
                numFields=numLeftwardFields;
            end
            
            for fii=1:numFields
                currFieldThetaData=currUnitStruct.thetaBasedTimeVsPhaseInfo.(currFieldDirStr).(sprintf('field%d',fii));

                allSpikesGoodLapIdxes=currFieldThetaData.allSpikeGoodLapIdxes;
                
                if(length(allSpikesGoodLapIdxes)<minNumSpikes)
                    continue
                end
                
                currFieldTimeBound=currUnitStruct.cycleTimeFieldBoundPerDirPerField{fii,di}-1; %-1 becaues manually entered before changing first theta cycle from 1 to 0
                
                currFieldTimeBound=currFieldTimeBound*0.9*0.9*0.95;
                
                allSpikesCurrTimeFieldFracs=currFieldThetaData.allLapSpikeCycleIDsFromFieldEntry/currFieldTimeBound;
                allSpikesCurrTimeFieldPhases=currFieldThetaData.allLapInFieldSpikePhases;
                
                goodLapSpikesCurrTimeFieldFracs=allSpikesCurrTimeFieldFracs(allSpikesGoodLapIdxes);
                goodLapSpikesCurrTimeFieldPhases=allSpikesCurrTimeFieldPhases(allSpikesGoodLapIdxes);
                
                if(maximizeInitializePhase)
                    %[bestShiftDeg,allPhasesShifted]=getBestShiftDeg1D(goodLapSpikesCurrTimeFieldPhases,goodLapSpikesCurrTimeFieldFracs,1/numTimeBins);
                    [bestShiftDeg,allPhasesShifted]=getBestShiftDeg1D(goodLapSpikesCurrTimeFieldPhases,goodLapSpikesCurrTimeFieldFracs,0.2);
                    
                    if(bestShiftDeg<30)
                        allPhases=allPhasesShifted;
                    else
                        %continue
                        allPhases=allPhasesShifted;
                    end
                else
                    allPhases=goodLapSpikesCurrTimeFieldPhases(:);
                end
               
            
                %{
                plot(goodLapSpikesCurrTimeFieldFracs,allPhases,'k.')
                hold on
                xlim([0 1])
                ylim([0 360])
                %}
                
                if(plotIndFields)
                    subplot(12,12,totalFieldCount+1)
                    plot(goodLapSpikesCurrTimeFieldFracs,allPhases,'k.','MarkerSize',1)
                    xlabel('Time in field (in cycles, frac.)')
                    ylabel('Theta phase (deg)')
                    xlim([0 1])
                    ylim([0 360])
                    
                end
                
                allUnitSpikeTimeFieldFracs=[allUnitSpikeTimeFieldFracs; goodLapSpikesCurrTimeFieldFracs(:)];
                allUnitSpikeTimeFieldPhases=[allUnitSpikeTimeFieldPhases; allPhases(:)];
                
                totalFieldCount=totalFieldCount+1;
            
            end
        end
        
    end
    
end
%%
fH=figure;

%{
numTimeBins=20;
numTimeBins=40;
numTimeBins=60;
numTimeBins=120;
numTimeBins=20;
%}

timeEdges=linspace(0,1,numTimeBins+1);
timeBinCenters=edgesToBins(timeEdges);
timeBinWidth=median(diff(timeBinCenters));

numPhaseBins=40;
numPhaseBins=60;
numPhaseBins=120;
%numPhaseBins=20;
phaseEdges=linspace(0,360,numPhaseBins+1);
phaseBinCenters=edgesToBins(phaseEdges);

[pxySmooth]=getJointDistr(allUnitSpikeTimeFieldFracs,allUnitSpikeTimeFieldPhases,timeEdges,phaseEdges,fH);


title(sprintf('n=%d fields',totalFieldCount))
xlabel('Time in field (in cycles, frac of field)')
                    ylabel('Theta phase (deg)')


daspect([1 360 1])
caxis([0 prctile(pxySmooth(:),97.5)])
%maxFig
setFigFontTo(18)
saveas(gcf,'timeInFieldVsThetaPhaseJointDistr.tif')

figure;
timeBinColors=jet(numTimeBins)
maxProbPhasePerTimeBin=NaN(numTimeBins,1);
for ti=1:numTimeBins
    currPdist=circSmooth(pxySmooth(ti,:),5);
    plot(phaseBinCenters,currPdist,'Color',timeBinColors(ti,:),'LineWidth',3)
    
    [~,maxID]=max(currPdist);
    maxProbPhasePerTimeBin(ti)=phaseBinCenters(maxID);
    hold on
end
colormap(gca,jet)
cb=colorbar
ylabel(cb,'Time in field (frac)')
xlabel('Theta phase (deg)')
ylabel('Probability')
xlim([0 360])

%maxFig
setFigFontTo(18)
saveas(gcf,'timeInFieldVsThetaPhaseDistr.tif')

figure;

plot(timeBinCenters,smooth(maxProbPhasePerTimeBin),'ko')
xlabel('Time in field (frac)')
ylabel('Most probable theta phase (deg)')
xlim([0 1])
ylim([0 360])

title(sprintf('n=%d fields',totalFieldCount))

daspect([1 360 1])
%maxFig
setFigFontTo(18)
saveas(gcf,'timeInFieldVsMaxProbThetaPhase.tif')

%autoArrangeFigures()
%{
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get circ mean phase per time bin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numTimeBinsForAvg=30;
timeAvgEdges=linspace(0,1,numTimeBinsForAvg+1);
timeAvgBinCenters=edgesToBins(timeAvgEdges);
timeAvgBinWidth=median(diff(timeAvgBinCenters));

phaseMeanPerTimeBin=NaN(numTimeBinsForAvg,1);

for ti=1:numTimeBinsForAvg
    currTimeBinStart=timeAvgEdges(ti)-timeAvgBinWidth/2;
    currTimeBinEnd=timeAvgEdges(ti)+timeAvgBinWidth/2;
    %currTimeBinCenter=timeBinCenters(ti);
    
    currBinIdxes=allUnitSpikeTimeFieldFracs>currTimeBinStart & allUnitSpikeTimeFieldFracs<=currTimeBinEnd;
    
    currBinPhases=allUnitSpikeTimeFieldPhases(currBinIdxes);
    currBinMeanPhase=circMeanDeg(currBinPhases);
    %currBinMeanPhase=circMedianDeg(currBinPhases);
    
    phaseMeanPerTimeBin(ti)=currBinMeanPhase;
    
end

figure;
plot(timeAvgBinCenters,phaseMeanPerTimeBin,'ko')
%}

close all; clear all; clc
dataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles';
cycleDataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData';
filePaths=getFilePathsRegex(dataDir,'*mat');
lowFreq=6;
highFreq=12;
numChannels=128;
chColors=jet(numChannels);
chColors=jet(10);
numPhaseBins=72;
phaseEdges=linspace(0,360,numPhaseBins+1);
phaseBinCenters=edgesToBins(phaseEdges);
startFileNum=1;

minSpeedForLocking=0.05;
%startFileNum=75;

showPlots=0;
showPlots=1;

for fi=startFileNum:length(filePaths)
       fi
    currFilePath=filePaths{fi};
    currFileName=getFileNameFromPath(currFilePath);
    fileBaseName=currFileName(1:(end-4));
    
    [currSessionNum,currSessName]=getSessionNum(fileBaseName);

    data=load(currFilePath);
    unitSpikeTimes=data.unitInfo.spikeTimes;
    thisUnitCh=data.unitInfo.thisUnitCh+10*(data.unitInfo.unitShankNum-1);
    
    spikePhasesPerCh=NaN(length(unitSpikeTimes),numChannels);
    phaseCountPerCh=NaN(numPhaseBins,numChannels);
    chKappas=NaN(numChannels,1);
    
    if(showPlots)
        figure; 
    end
    
    for ci=1:numChannels
        if(ci<65)
            sideNum=1;
            titleStr='Left hemisphere CA1 LFP';
        else
             sideNum=2;
             titleStr='Right hemisphere CA1 LFP';
        end
            cycleInfoFileName=fullfile(cycleDataDir,sprintf('%s_Ch%d_%d-%d_minmax.mat',currSessName,ci,lowFreq,highFreq));
            currUnitCycleInfo=load(cycleInfoFileName);
            
            minTimes=currUnitCycleInfo.cycleMinTimes;
                maxTimes=currUnitCycleInfo.cycleMaxTimes;
                
                absSpeedPerTime=abs(data.filledSignedSpeedPerTime);
                positionTimeAxis=data.positionTimeAxis;
                
              unitSpikeSpeeds=interp1(positionTimeAxis,absSpeedPerTime,unitSpikeTimes);
                
            
            numCycles=length(currUnitCycleInfo.cycleMaxTimes);
                [spikeAssignMax spikeAssignMin spikeAssignZero] = assignSpikesToCycles2017(unitSpikeTimes, minTimes, maxTimes);

                mode=2; %phase from max to max
                %mode=1;%phase from startMin to endMin
                calcTimeOffsets=1;
                [currChSpikePhases spikeOffsets spikeOffsetsEnd] = assignPhaseToSpikes(unitSpikeTimes, spikeAssignMax, spikeAssignMin, minTimes, maxTimes, mode, calcTimeOffsets, spikeAssignZero);
            spikePhasesPerCh(:,ci)=currChSpikePhases;
            
            
            nonNanSpikePhases=currChSpikePhases(~isnan(currChSpikePhases) & (unitSpikeSpeeds(:)')>minSpeedForLocking);
            chKappas(ci)=circ_kappa(deg2rad(nonNanSpikePhases));
            
            %phaseCountPerCh(:,ci)=histcounts(currChSpikePhases,phaseEdges);
            
            phaseCountPerCh(:,ci)=histcounts(nonNanSpikePhases,phaseEdges);
            if(showPlots)
                colorNum=mod(ci-1,10)+1;
                subplot(2,2,sideNum)
                plot(phaseBinCenters,smooth(phaseCountPerCh(:,ci)),'Color',chColors(colorNum,:),'LineWidth',2)
                title(titleStr)
                hold on
                xlim([0 360])
                xlabel('Spike phase')
                ylabel('Count')
                colormap(jet)
                cb=colorbar;
                ylabel(cb,'Channel num mod 10')
                caxis([0 9])
                hold on
            end
            %drawnow
            ci
    end
    
    if(showPlots)
        subplot(2,2,[3 4])

       plot(chKappas)
       hold on
        plot(chKappas,'k.','MarkerSize',20)
       plot([64.5 64.5],ylim,'k--','LineWidth',5)
       chH=plot([thisUnitCh thisUnitCh],ylim,'r-','LineWidth',5)
       legend(chH,'cluster ch. num','Location','best')
       xlim([1 numChannels])
       xlabel('Ch num')
       ylabel('Phase locking (kappa)')

        uberTitle(removeUnderscores(sprintf('%s theta phase across %d channels',fileBaseName,numChannels)))
        maxFig
        setFigFontTo(18)
        saveas(gcf,sprintf('%s_SelectMasterTheta.tif',fileBaseName))

    end
    
    chKappasZscore=zscoreLFP(chKappas);
    save( currFilePath, 'spikePhasesPerCh','phaseCountPerCh','chKappas','chKappasZscore','-append')
    if(showPlots)
        close all
    end
  end       
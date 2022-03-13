close all; clear all; clc
dataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData';
saveDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/placeCellQualityImages';

%{
Grosberg & Buzsaki, 2016 place field analysis method:
For each well isolated principal cell a spike firing-by-position vector was
constructed by 
1) binning its spikes in non-overlapping 2 cm bins. 
2) this vector was smoothed with a 5 cm Gaussian kernel, and
3) divided by the smoothed (5 cm Gaussian kernel) occupancy-by-position 
resulting in a smoothed position by-firing rate vector. 

In
the case of the circular maze location was linearized and defined as starting at the edge of
reward area, and increasing clockwise, terminating at the opposing edge of the reward
area. 
The hypothesis of place-selective firing was tested by constructing 5,000 null firing
rate vectors in which the cell?s spikes were randomly sampled using the un-smoothed
occupancy-by-position vector as the probability density function, smoothed with a 5cm
Gaussian kernel and divided by the smoothed occupancy-by-position vector. 

A cell was determined to have a place field (and thus, to be a place cell) if at least 5 consecutive bins
were above the 99th percentile of their null distributions and the cell exhibited a withinfield peak firing rate of at least 1Hz. All place field detection analysis was restricted to
epochs during which the animal?s velocity was at least 5 cm/s and in which the rat was
outside of the reward areas. 

For linear track maze runs all place field analysis was
carried out independently for left and right directions of movement (3, 4, 18). Standard
measures of place field properties (information per spike, information per second, spatial
coherence and stability) were performed in the same manner as in (32)). 
%}

filePaths=getFilePathsRegex(dataDir,'*Info.mat');

for fi=1:length(filePaths)
    tic
    currFilePath=filePaths{fi};
    
    fileName=getFileNameFromPath(currFilePath);
    baseName=fileName(1:(end-4));
    
    if(~strContainsCircSessionName(fileName))
        continue
    end
   
    
    
    currPlaceInfo=load(currFilePath);
    currUnitInfo=currPlaceInfo.unitInfo;
     baseName=[baseName '_'  currUnitInfo.unitCellType];
     
    mazeTaskTimeBounds=currUnitInfo.mazeTaskTimeBounds;
    
    firingRateOverTime=currUnitInfo.firingRateOverTime;
    rateTimeBinCenters=currUnitInfo.rateTimeBinCenters;
    spikePhases=currUnitInfo.spikePhases;
    spikeTimes=currUnitInfo.spikeTimes;
    
   positionTimeAxis=currPlaceInfo.positionTimeAxis;
   positionPerTimeStep=currPlaceInfo.positionPerTimeStep;
   spikePositions=interp1(positionTimeAxis,positionPerTimeStep,spikeTimes);
    
    
    [positionBins,firingRatePerPositionRight, firingRatePerPositionLeft,spikeTimeDirectionAssignment] = getSpikeStatsPerPosition(spikeTimes,positionTimeAxis,positionPerTimeStep,currPlaceInfo.speedPerTimeStep);
   
    peakRateRight=max(firingRatePerPositionRight);
    peakRateLeft=max(firingRatePerPositionLeft);
    
    
    rightSpikeTimes=spikeTimes(spikeTimeDirectionAssignment==1);
     leftSpikeTimes=spikeTimes(spikeTimeDirectionAssignment==2);
     
     rightSpikePhases=spikePhases(spikeTimeDirectionAssignment==1);
     leftSpikePhases=spikePhases(spikeTimeDirectionAssignment==2);
     
     rightSpikePositions=spikePositions(spikeTimeDirectionAssignment==1);
     leftSpikePositions=spikePositions(spikeTimeDirectionAssignment==2);
    
    if(max([peakRateRight peakRateLeft])>5)
        figure;
        subplot(2,2,1)
        plot(positionBins,firingRatePerPositionRight,'k-','LineWidth',4)
               title('right lap')
               xlim([min(positionBins) max(positionBins)])
                xlabel('Position in track (m)')
                ylabel('Firing rate (Hz)')
                 box off
        subplot(2,2,2)
        hold on
          plot(positionBins,firingRatePerPositionLeft,'k-','LineWidth',4)
   xlim([min(positionBins) max(positionBins)])
    xlabel('Position in track (m)')
   ylabel('Firing rate (Hz)')
          title('left lap')
          
           box off
         subplot(2,2,3)
        plot(rightSpikePositions,rightSpikePhases,'k.')
        xlim([min(positionBins) max(positionBins)])
           ylim([0 360])
            ylabel('Theta phase (deg)')
             xlabel('Position in track (m)')
             box off
        subplot(2,2,4)
        hold on
        plot(leftSpikePositions,leftSpikePhases,'k.')
        xlim([min(positionBins) max(positionBins)])
        ylim([0 360])
        ylabel('Theta phase (deg)')
        
        xlabel('Position in track (m)')
         box off
        
        
        
        uberTitle(removeUnderscores(baseName))
        maxFig
        setFigFontTo(18)
         saveas(gcf,fullfile(saveDir,sprintf('%s_PlaceFieldAndPhaseVsPos.tif',baseName)))
           %disp('')
           close all
    end
    
 
    save(currFilePath,'positionBins','firingRatePerPositionRight','firingRatePerPositionLeft','spikeTimeDirectionAssignment','peakRateRight','peakRateLeft',...
        'rightSpikeTimes','leftSpikeTimes','rightSpikePhases','leftSpikePhases','leftSpikePositions','rightSpikePositions' ,'-append');
    
    
    data=load(currFilePath)
    toc
end

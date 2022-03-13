function [positionBinCenters,firingRatePerPositionRight, firingRatePerPositionLeft,...
	timeInLapBins,firingRatePerTimeRight, firingRatePerTimeLeft]...
	 = getFiringRatePerPosition(spikeTimes,positionTimeAxis,positionPerTimeStep,timeInLapPerTimeStep,speedPerTimeStep,minPos,maxPos)
%{
Grosberg & Buzsaki, 2016 place field analysis method:
For each well isolated principal cell a spike firing-by-position vector was
constructed by 
1) binning its spikes in non-overlapping 2 cm bins. 
2) this vector was smoothed with a 5 cm Gaussian kernel, and
3) divided by the smoothed (5 cm Gaussian kernel) occupancy-by-position 
resulting in a smoothed position by-firing rate vector. 

A cell was determined to have a place field (and thus, to be a place cell) if at least 5 consecutive bins
were above the 99th percentile of their null distributions and the cell exhibited a withinfield peak firing rate of at least 1Hz. 
All place field detection analysis was restricted to
epochs during which the animal?s velocity was at least 5 cm/s and in which the rat was
outside of the reward areas. 
%}

positionPerTimeStepRight=positionPerTimeStep;
%positionPerTimeStepRight(speedPerTimeStep<0)=NaN;
positionPerTimeStepRight(speedPerTimeStep<0.05)=NaN;

positionPerTimeStepLeft=positionPerTimeStep;
%positionPerTimeStepLeft(speedPerTimeStep>0)=NaN;
positionPerTimeStepLeft(speedPerTimeStep>-0.05)=NaN;

timeInLapPerTimeStepRight=timeInLapPerTimeStep;
%timeInLapPerTimeStepRight(speedPerTimeStep<0)=NaN;
timeInLapPerTimeStepRight(speedPerTimeStep<0.05)=NaN;

timeInLapPerTimeStepLeft=timeInLapPerTimeStep;
%timeInLapPerTimeStepLeft(speedPerTimeStep>0)=NaN;
timeInLapPerTimeStepLeft(speedPerTimeStep>-0.05)=NaN;

approxTimeStep=median(diff(positionTimeAxis));

%positionPerTimeStepRight=fillmissing(positionPerTimeStepRight,'linear');
%positionPerTimeStepLeft=fillmissing(positionPerTimeStepLeft,'linear');

%{
figure; 
plot(positionTimeAxis,positionPerTimeStepRight,'k','LineWidth',3)
figure
plot(positionTimeAxis,positionPerTimeStepLeft,'b','LineWidth',3)
%}

posBinWidth=0.02;
posBinWidth=0.03;
positionBinEdges=min(positionPerTimeStep):posBinWidth:max(positionPerTimeStep);
positionBinCenters=edgesToBins(positionBinEdges);

outOfFieldPosBins=positionBinCenters<minPos | positionBinCenters>maxPos;

timeBinWidth=0.2;
maxTimeInLap=2.5;

timeBinWidth=0.1;
maxTimeInLap=3;

timeBinWidth=0.05;
maxTimeInLap=5;
maxTimeInLap=5;

timeBinWidth=0.075;
%timeBinWidth=0.1;
%timeBinWidth=0.0625;

timeInLapBinEdges=0:timeBinWidth:maxTimeInLap;
timeInLapBins=edgesToBins(timeInLapBinEdges);

currDirSpikeCountPerPos=zeros(length(positionBinCenters),1);
currDirSpikeCountPerTimeInLap=zeros(length(timeInLapBins),1);

currDirTimeCountPerTimeInLap=zeros(length(timeInLapBins),1);
currDirTimeCountPerPos=zeros(length(positionBinCenters),1);
for direction=1:2
    if(direction==1)
        currPositionPerTime=positionPerTimeStepRight;
        currTimeInLapPerTime=timeInLapPerTimeStepRight;
    else
	currPositionPerTime=positionPerTimeStepLeft;
        currTimeInLapPerTime=timeInLapPerTimeStepLeft;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %1) binning its spikes in non-overlapping 2 cm bins
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    positionPerTimeAsBinNums=discretize(currPositionPerTime,positionBinEdges);
    timeInLapPerTimeAsBinNums=discretize(currTimeInLapPerTime,timeInLapBinEdges);
    
    %count spikes corresponding to each position bin
    for ti=1:length(positionTimeAxis)
        
        currPosBin=positionPerTimeAsBinNums(ti);
        currTimeInLapBin=timeInLapPerTimeAsBinNums(ti);
        
        if(currPositionPerTime(ti)<minPos || currPositionPerTime(ti)>maxPos)
            continue
        end 
        
        if(isnan(currPosBin))
            continue
        end
        
         if(isnan(currTimeInLapBin))
            continue
        end
        
        binStartTime=positionTimeAxis(ti)-approxTimeStep/2;
        binEndTime=positionTimeAxis(ti)+approxTimeStep/2;
        
        spikeTimesInBin=find(spikeTimes>=binStartTime & spikeTimes<binEndTime);
        
        currDirSpikeCountPerPos(currPosBin)=currDirSpikeCountPerPos(currPosBin)+length(spikeTimesInBin);
        currDirSpikeCountPerTimeInLap(currTimeInLapBin)=currDirSpikeCountPerTimeInLap(currTimeInLapBin)+length(spikeTimesInBin);
      	currDirTimeCountPerTimeInLap(currTimeInLapBin)=currDirTimeCountPerTimeInLap(currTimeInLapBin)+approxTimeStep;  
	    currDirTimeCountPerPos(currPosBin)=currDirTimeCountPerPos(currPosBin)+approxTimeStep;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %2) this vector was smoothed with a 5 cm Gaussian kernel
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    smoothBinWidthPos=3;
    smoothBinWidthPos=5;
    smoothBinWidthTime=5;
     smoothBinWidthPos=7;
    smoothBinWidthTime=7;
     smoothBinWidthPos=8;
    smoothBinWidthTime=8;
    %smoothBinWidthPos=6;
    %smoothBinWidthTime=6;
    
    currDirSpikeCountPerPos=smoothdata(currDirSpikeCountPerPos(:),'gaussian',smoothBinWidthPos); %3 bin wide window (6cm)
    currDirTimeCountPerPos=smoothdata(currDirTimeCountPerPos(:),'gaussian',smoothBinWidthPos);
    
    currDirSpikeCountPerTimeInLap=smoothdata(currDirSpikeCountPerTimeInLap(:),'gaussian',smoothBinWidthTime); %sec wide window
    currDirTimeCountPerTimeInLap=smoothdata(currDirTimeCountPerTimeInLap(:),'gaussian',smoothBinWidthTime);
    
    %{
    figure;
    
    plot(positionBins,currDirSpikeCountPerPos,'k-','LineWidth',3)
    hold on
    plot(positionBins,currDirTimeCountPerPos,'b-','LineWidth',3)
    %}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %3) divide by the smoothed (5 cm Gaussian kernel) occupancy-by-position 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    if(direction==1)
        firingRatePerPositionRight=currDirSpikeCountPerPos./currDirTimeCountPerPos;
   	    firingRatePerTimeRight=currDirSpikeCountPerTimeInLap./currDirTimeCountPerTimeInLap; 
        
        firingRatePerPositionRight(outOfFieldPosBins)=NaN;
        %firingRatePerPositionRight(outOfFieldPosBins)=0;
         %yyaxis right
         %plot(positionBins,firingRatePerPositionRight,'r-','LineWidth',3)
    else
        firingRatePerPositionLeft=currDirSpikeCountPerPos./currDirTimeCountPerPos;
   	    firingRatePerTimeLeft=currDirSpikeCountPerTimeInLap./currDirTimeCountPerTimeInLap; 
        firingRatePerPositionLeft(outOfFieldPosBins)=NaN;
        %firingRatePerPositionLeft(outOfFieldPosBins)=0;
        %yyaxis right
         %plot(positionBins,firingRatePerPositionLeft,'r-','LineWidth',3)
    end
   
   
end



end


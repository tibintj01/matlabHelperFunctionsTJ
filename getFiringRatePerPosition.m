function [positionBins,firingRatePerPositionRight, firingRatePerPositionLeft] = getFiringRatePerPosition(spikeTimes,positionTimeAxis,positionPerTimeStep,speedPerTimeStep)
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

approxTimeStep=median(diff(positionTimeAxis));

%positionPerTimeStepRight=fillmissing(positionPerTimeStepRight,'linear');
%positionPerTimeStepLeft=fillmissing(positionPerTimeStepLeft,'linear');

%{
figure; 
plot(positionTimeAxis,positionPerTimeStepRight,'k','LineWidth',3)
figure
plot(positionTimeAxis,positionPerTimeStepLeft,'b','LineWidth',3)
%}

positionBinEdges=min(positionPerTimeStep):0.02:max(positionPerTimeStep);
positionBins=edgesToBins(positionBinEdges);

currDirSpikeCountPerPos=zeros(length(positionBins),1);
currDirTimeCountPerPos=zeros(length(positionBins),1);
for direction=1:2
    if(direction==1)
        currPositionPerTime=positionPerTimeStepRight;
    else
        currPositionPerTime=positionPerTimeStepLeft;
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %1) binning its spikes in non-overlapping 2 cm bins
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    positionPerTimeAsBinNums=discretize(currPositionPerTime,positionBinEdges);
    
    %count spikes corresponding to each position bin
    for ti=1:length(positionTimeAxis)
        
        currPosBin=positionPerTimeAsBinNums(ti);
        
        if(isnan(currPosBin))
            continue
        end
        
        binStartTime=positionTimeAxis(ti)-approxTimeStep/2;
        binEndTime=positionTimeAxis(ti)+approxTimeStep/2;
        
        spikeTimesInBin=find(spikeTimes>=binStartTime & spikeTimes<binEndTime);
        
        currDirSpikeCountPerPos(currPosBin)=currDirSpikeCountPerPos(currPosBin)+length(spikeTimesInBin);
        currDirTimeCountPerPos(currPosBin)=currDirTimeCountPerPos(currPosBin)+approxTimeStep;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %2) this vector was smoothed with a 5 cm Gaussian kernel
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    currDirSpikeCountPerPos=smoothdata(currDirSpikeCountPerPos(:),'gaussian',3); %3 bin wide window (6cm)
    currDirTimeCountPerPos=smoothdata(currDirTimeCountPerPos(:),'gaussian',3);
    
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
         %yyaxis right
         %plot(positionBins,firingRatePerPositionRight,'r-','LineWidth',3)
    else
        firingRatePerPositionLeft=currDirSpikeCountPerPos./currDirTimeCountPerPos;
        %yyaxis right
         %plot(positionBins,firingRatePerPositionLeft,'r-','LineWidth',3)
    end
   
    
end



end


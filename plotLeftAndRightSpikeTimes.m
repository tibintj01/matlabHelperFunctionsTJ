function [] = plotLeftAndRightSpikeTimes(unitData)


rightwardSpikeTimes=unitData.rightSpikeTimes;
leftwardSpikeTimes=unitData.leftSpikeTimes;


rightwardSpikePositions=unitData.rightSpikePositions;
leftwardSpikePositions=unitData.leftSpikePositions;

if(~isempty(rightwardSpikePositions))
    plot(rightwardSpikeTimes, rightwardSpikePositions,'m*','MarkerSize',10) 
    hold on
end
if(~isempty(leftwardSpikePositions))
    hold on
plot(leftwardSpikeTimes,leftwardSpikePositions,'b*','MarkerSize',10) 
end
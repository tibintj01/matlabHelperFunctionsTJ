function [spikeLaps]=getLapPerSpikeTime(spikeTimes,lapStartTimes,lapEndTimes)
%return lap of each spike given spike times and lap start and end times
spikeLaps=NaN(size(spikeTimes));
numSpikes=length(spikeTimes);

numLaps=length(lapStartTimes);

for si=1:numSpikes
    currSpikeTime=spikeTimes(si);
    currSpikeIsInLap=currSpikeTime>=lapStartTimes & currSpikeTime<lapEndTimes;
    
    currSpikeLap=find(currSpikeIsInLap);
    if(~isempty(currSpikeLap))
        spikeLaps(si)=currSpikeLap;
    end
end



figure; 

inFieldCyclesPerSpikeGood=inFieldCyclesPerSpike(~inBadLap);
inFieldDistsPerSpikeGood=inFieldDistsPerSpike(~inBadLap);

numPts=length(inFieldCyclesPerSpikeGood);

for dispIdx=numPts-1:numPts
    plot(inFieldCyclesPerSpikeGood(1:dispIdx), inFieldDistsPerSpikeGood(1:dispIdx),'k')
    hold on
    plot(inFieldCyclesPerSpikeGood(1:dispIdx), inFieldDistsPerSpikeGood(1:dispIdx),'ko')
    
    drawnow
    pause(.1)
    %close all
    
end
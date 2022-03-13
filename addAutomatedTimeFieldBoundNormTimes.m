data=load('allUnitInFieldGoodSpeedSpikePhasesAndTimes.mat')


numFields=length(data.spikePhasesPerField);

autoRenormedSpikeTimeFracInFieldPerField=cell(numFields,1);

for fi=1:numFields
    
    
    currFieldPhases=data.spikePhasesPerField{fi};
    currFieldTimeFracs=data.spikeTimeFracInFieldPerField{fi};
    
    if(isempty(currFieldPhases) || isempty(currFieldTimeFracs))
        continue
    end
    
    currFieldRenormedSpikeTimeFracs=getRealTimeBoundNormedSpikeFracs(currFieldPhases,currFieldTimeFracs);
    
    autoRenormedSpikeTimeFracInFieldPerField{fi}=currFieldRenormedSpikeTimeFracs;
end

save('allUnitInFieldGoodSpeedSpikePhasesAndTimes.mat','autoRenormedSpikeTimeFracInFieldPerField','-append')

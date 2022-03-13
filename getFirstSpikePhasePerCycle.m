function [spikesPerCycleInfo]=getFirstSpikePhasePerCycle(cycleInfo,spikePhase,spikeAssignMax)
%UNTITLED5 Summary of this function goes here
    
    %cycleNums=unique(spikeAssignMax)
    %cycleMaxTimes=cycleInfo.cycleMaxTimes;
    numCycles=length(cycleInfo.cycleMaxTimes);
   
    allSpikesPerCycle=NaN(numCycles,150); 
    firstSpikePerCycle=NaN(numCycles,1);
    linearAdaptationWeightedCOMPhase=NaN(numCycles,1);
    expAdaptationWeightedCOMPhase=NaN(numCycles,1);
    numSpikesPerCycle=zeros(numCycles,1);
    for currCycleNum=1:numCycles
        %currCycleNum=cycleNums(currCycleNumID);
        spikeIDsInCurrCycle=find(currCycleNum==spikeAssignMax);
        spikePhasesInCurrCycle=spikePhase(spikeIDsInCurrCycle);

      
        if(~isempty(spikePhasesInCurrCycle))
	    allSpikesPerCycle(currCycleNum,1:length(spikePhasesInCurrCycle))=spikePhasesInCurrCycle;  
            firstSpikePerCycle(currCycleNum)=spikePhasesInCurrCycle(1);
            numSpikesPerCycle(currCycleNum)=length(spikePhasesInCurrCycle);
            linearMassWeightedTotal=0;
            expMassWeightedTotal=0;
            for i=1:length(spikePhasesInCurrCycle)
                currLinearWeight=length(spikePhasesInCurrCycle)-i+1;
                currExpWeight=exp(currLinearWeight);
                
                linearMassWeightedTotal=linearMassWeightedTotal+currLinearWeight*spikePhasesInCurrCycle(i);
                expMassWeightedTotal=expMassWeightedTotal+currExpWeight*spikePhasesInCurrCycle(i);

            end
            totalLinearMass=sum(1:length(spikePhasesInCurrCycle));
            totalExpMass=sum(exp(1:length(spikePhasesInCurrCycle)));
            
            linearAdaptationWeightedCOMPhase(currCycleNum)=linearMassWeightedTotal/totalLinearMass;
            expAdaptationWeightedCOMPhase(currCycleNum)=expMassWeightedTotal/totalExpMass;
        end
    end
    
    spikesPerCycleInfo.cycleStartTimes=cycleInfo.cycleMaxTimes;
    spikesPerCycleInfo.firstSpikePhasePerCycle=firstSpikePerCycle(:);
    spikesPerCycleInfo.linearAdaptationWeightedCOMPhase=linearAdaptationWeightedCOMPhase(:);
    spikesPerCycleInfo.expAdaptationWeightedCOMPhase=expAdaptationWeightedCOMPhase(:);
    spikesPerCycleInfo.numSpikesPerCycle=numSpikesPerCycle(:);
    spikesPerCycleInfo.allSpikesPerCycle=allSpikesPerCycle;
end


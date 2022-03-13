if(di==2)
            currFieldDirection='leftward';
            goingBackwardIdxes=data.speedPerTimeStep>0;
            gaussSpikeRatePerTime=data.leftLocalGaussSpikeRatePerPosTime;
            refUnitData=leftwardRefUnitData;

            if(useRefThetaCh && ~isnan(data.refThetaCh))
             spikePhases=data.leftSpikesRefPhasesMaxToMax;
            else
             spikePhases=data.leftSpikePhases;

            end

            spikeTimes=data.leftSpikeTimes;
            spikePositions=data.leftSpikePositions;

            if(justFirstSpikes)
                if(~isfield(data,'leftIsFirstSpikeInRefCycle'))
                    continue
                end
                spikeTimes=spikeTimes(data.leftIsFirstSpikeInRefCycle);
                spikePositions=spikePositions(data.leftIsFirstSpikeInRefCycle);
                spikePhases=spikePhases(data.leftIsFirstSpikeInRefCycle);
            end

        elseif(di==1)
            if(strContainsCircSessionName(currFileName))
                continue %circular mazes always leftward laps
            end
            currFieldDirection='rightward';
            goingBackwardIdxes=data.speedPerTimeStep<0;
            refUnitData=rightwardRefUnitData;

            gaussSpikeRatePerTime=data.rightLocalGaussSpikeRatePerPosTime;

            if(useRefThetaCh && ~isnan(data.refThetaCh))
              spikePhases=data.rightSpikesRefPhasesMaxToMax;
            else
             spikePhases=data.rightSpikePhases;

            end

            spikeTimes=data.rightSpikeTimes;
            spikePositions=data.rightSpikePositions;

            if(justFirstSpikes)
                if(~isfield(data,'rightIsFirstSpikeInRefCycle'))
                    continue
                end
                spikeTimes=spikeTimes(data.rightIsFirstSpikeInRefCycle);
                spikePositions=spikePositions(data.rightIsFirstSpikeInRefCycle);
                spikePhases=spikePhases(data.rightIsFirstSpikeInRefCycle);
            end
        end

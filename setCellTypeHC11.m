cellType=NaN;
        if(ismember(currIDnum,spikeInfo.sessInfo.Spikes.PyrIDs))
            cellType='pyr';
        elseif(ismember(currIDnum,spikeInfo.sessInfo.Spikes.IntIDs))
            cellType='int';
        end
 unitThetaBasedDataStruct=dataStruct.thetaBasedTimeVsPhaseInfo.(currFieldDirStr);
        currFieldSpikeTimesInExp=spikeTimesInExpPerField{fi};
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %recover field number
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        containsSpikeTol=0.001;
        numFieldsInStruct=length(fieldnames(unitThetaBasedDataStruct));
        
        if(numFieldsInStruct>1)
            
            disp('')
            
        end
        for f=1:numFieldsInStruct
            spikeTimesThisFieldNumCandidate=unitThetaBasedDataStruct.(sprintf('field%d',f)).allLapInFieldSpikeTimesInExp;
            
            allSpikesContained=1;
            %if theta based field does not contain one of these spikes it
            %is not the same field
            for s=1:length(currFieldSpikeTimesInExp)
                currSpikeTime=currFieldSpikeTimesInExp(s);
                if(min(abs(currSpikeTime-spikeTimesThisFieldNumCandidate))>containsSpikeTol)
                    allSpikesContained=0;
                end
            end
            
            if(allSpikesContained==1)
                currFieldThetaData=unitThetaBasedDataStruct.(sprintf('field%d',f));
                fieldNumWithinUnit=f;
                break
            end
        
        end
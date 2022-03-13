
representativeLeftwardUnitInfo=[];
representativeRightwardUnitInfo=[];

    for ri=1:length(currSessUnitFilePaths)
        representativeUnitInfoCandidate=load(currSessUnitFilePaths{ri});
        
        if(~isempty(representativeUnitInfoCandidate.rightwardFieldStartEndM))
            representativeRightwardUnitInfo=representativeUnitInfoCandidate;
        end
        
        if(~isempty(representativeUnitInfoCandidate.leftwardFieldStartEndM))
            representativeLeftwardUnitInfo=representativeUnitInfoCandidate;
        end
        
        if(~isempty(representativeLeftwardUnitInfo) && ~isempty(representativeRightwardUnitInfo))
            break
        end
        
    end
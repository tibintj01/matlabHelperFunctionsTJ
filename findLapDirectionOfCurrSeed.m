  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %shift back start pointer until its beyond touchdown area
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        firstIdxOfThisLap=currLapSeedPosIdx;
        while(abs(filledPositionPerTime(firstIdxOfThisLap)-maxPositionHigh)>closenessThresh && abs(filledPositionPerTime(firstIdxOfThisLap)-maxPositionLow)>closenessThresh)
            if(firstIdxOfThisLap==1)
               break
           end
            firstIdxOfThisLap=firstIdxOfThisLap-1;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %shift fwd end pointer until its beyond touchdown area
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        lastIdxOfThisLap=currLapSeedPosIdx;
        while(abs(filledPositionPerTime(lastIdxOfThisLap)-maxPositionHigh)>closenessThresh && abs(filledPositionPerTime(lastIdxOfThisLap)-maxPositionLow)>closenessThresh)
           if(lastIdxOfThisLap==length(filledPositionPerTime))
               break
           end
            lastIdxOfThisLap=lastIdxOfThisLap+1;
        end
        
        if(filledPositionPerTime(lastIdxOfThisLap) < filledPositionPerTime(firstIdxOfThisLap) )
            currSeedInLeftwardLap=1;
        else
            currSeedInLeftwardLap=0;
        end
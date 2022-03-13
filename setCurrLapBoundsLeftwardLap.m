        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %shift back start pointer until its beyond touchdown area
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        firstIdxOfThisLap=currLapSeedPosIdx;
        while(filledPositionPerTime(firstIdxOfThisLap)<maxPositionHigh)
            if(firstIdxOfThisLap==1)
               break
           end
            firstIdxOfThisLap=firstIdxOfThisLap-1;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %shift back start pointer until speed crosses zero
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        while(filledSignedSpeedPerTime(firstIdxOfThisLap)<0)
            if(firstIdxOfThisLap==1)
               break
           end
            firstIdxOfThisLap=firstIdxOfThisLap-1;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %shift fwd start pointer until its beyond touchdown area
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        lastIdxOfThisLap=currLapSeedPosIdx;
        while(filledPositionPerTime(lastIdxOfThisLap)>maxPositionLow)
            if(lastIdxOfThisLap==length(filledPositionPerTime))
               break
           end
            lastIdxOfThisLap=lastIdxOfThisLap+1;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %shift fwd start pointer until speed crosses zero
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        while(filledSignedSpeedPerTime(lastIdxOfThisLap)<0)
            if(lastIdxOfThisLap==length(filledPositionPerTime))
               break
           end
            lastIdxOfThisLap=lastIdxOfThisLap+1;
        end

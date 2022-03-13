function [startBoundTime,endBoundTime,withinBoundsTimeIdxes] = getTimeBoundsOutwardFromSeedAndPosBounds(lowestPossibleTime,highestPossibleTime,timeAxis,posPerTime,centralSeedTime,lowerPosBound,upperPosBound)
%UNTITLED Summary of this function goes here

         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %step backward from current cycle seed to find start of trav
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         [~,lowestPossibleTimeIdx]=min(abs(lowestPossibleTime - timeAxis));
         
         [~,centralSeedTimeIdxToStart]=min(abs(centralSeedTime - timeAxis));
         
         while(posPerTime(centralSeedTimeIdxToStart)>lowerPosBound)
             centralSeedTimeIdxToStart=centralSeedTimeIdxToStart-1;
             
             if(centralSeedTimeIdxToStart<lowestPossibleTimeIdx)
                 break
             end
         end
         
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %step forward from current cycle seed to find end of trav
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         [~,highestPossibleTimeIdx]=min(abs(highestPossibleTime - timeAxis));
         [~,centralSeedTimeIdxToEnd]=min(abs(centralSeedTime - timeAxis));
         
         while(standardizedPosPerTime(centralSeedTimeIdxToEnd)<upperPosBound)
             centralSeedTimeIdxToEnd=centralSeedTimeIdxToEnd+1;
             
             if(centralSeedTimeIdxToEnd>highestPossibleTimeIdx)
                 break
             end
         end
         

         withinBoundsTimeIdxes=centralSeedTimeIdxToStart:centralSeedTimeIdxToEnd;
         
         startBoundTime=timeAxis(centralSeedTimeIdxToStart);
         endBoundTime=timeAxis(currTravPosTimeIdxEnd);
         

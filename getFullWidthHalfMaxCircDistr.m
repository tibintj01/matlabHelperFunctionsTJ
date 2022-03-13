function [currDistrWidth,risingIntersectPhase,risingIntersectProb,descendingIntersectPhase,descendingIntersectProb]=getFullWidthHalfMaxCircDistr(phaseBinCenters,currPdist)
%finds full width half max of circ distribution using angDiffDeg
%making sure rising is before descending so that peak is captured

    phaseBinDiff=median(diff(phaseBinCenters));

currPdistOriginal=currPdist(:);

upSampleFactor=10;
upSampleFactor=30;

[phaseBinCenters,currPdist,padLength]=getCircPaddedDistr(phaseBinCenters,currPdist,upSampleFactor);


    %[maxProb,maxID]=max(currPdistOriginal);
    %maxID=maxID+padLength;
    
    [maxProb,maxID]=max(currPdist((padLength+1):(end-padLength)));
    [minProb,minID]=min(currPdist((padLength+1):(end-padLength)));
    
    maxID=maxID+padLength;

    %halfMaxProb=minProb+(maxProb-minProb)/2;
     halfMaxProb=maxProb/2;
    
    
    risingIntersectPhase=NaN;
    descendingIntersectPhase=NaN;
    
    risingIntersectI=maxID;
    while(currPdist(risingIntersectI)>halfMaxProb)
       
        
        if(risingIntersectI==1)
            break
        end
         risingIntersectI=risingIntersectI-1;
    end
    
    risingIntersectPhase=phaseBinCenters(risingIntersectI);
    risingIntersectProb=currPdist(risingIntersectI);
    
   
    
    fallingIntersectI=maxID;
    while(currPdist(fallingIntersectI)>halfMaxProb)
       
        
        if(fallingIntersectI==length(currPdist))
            break
        end
         fallingIntersectI=fallingIntersectI+1;
    end
       descendingIntersectPhase=phaseBinCenters(fallingIntersectI);
      descendingIntersectProb=currPdist(fallingIntersectI);
      
      currDistrWidth=(fallingIntersectI-risingIntersectI)*phaseBinDiff/upSampleFactor;
      
      if(currDistrWidth==0 || currDistrWidth>360) %not interpretable, width detection failed
          currDistrWidth=NaN;
      end
      
    
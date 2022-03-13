function [currDistrWidth,risingIntersectPhase,risingIntersectProb,descendingIntersectPhase,descendingIntersectProb]=getFullWidthHalfMaxLinearDistr(xBinCenters,currPdist)
%finds full width half max of circ distribution using angDiffDeg
%making sure rising is before descending so that peak is captured

    xBinDiff=median(diff(xBinCenters));

    currPdistOriginal=currPdist(:);

    upSampleFactor=10;
    upSampleFactor=30;
      %upSampleFactor=100;
      
      upsampledXbinCenters=min(xBinCenters):(xBinDiff/upSampleFactor):max(xBinCenters);
      
      currPdist=interp1(xBinCenters,currPdist,upsampledXbinCenters);
    
      xBinCenters=upsampledXbinCenters;
      
    [maxProb,maxID]=max(currPdist);


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
    
    risingIntersectPhase=xBinCenters(risingIntersectI);
    risingIntersectProb=currPdist(risingIntersectI);
    
   
    
    fallingIntersectI=maxID;
    while(currPdist(fallingIntersectI)>halfMaxProb)
       
        
        if(fallingIntersectI==length(currPdist))
            break
        end
         fallingIntersectI=fallingIntersectI+1;
    end
       descendingIntersectPhase=xBinCenters(fallingIntersectI);
      descendingIntersectProb=currPdist(fallingIntersectI);
      
      currDistrWidth=(fallingIntersectI-risingIntersectI)*xBinDiff/upSampleFactor;
      
      if(currDistrWidth==0) %not interpretable, width detection failed
          currDistrWidth=NaN;
      end
      
    
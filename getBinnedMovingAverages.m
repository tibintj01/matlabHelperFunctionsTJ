function [yMovAvgValues,yMovSEMValues] = getBinnedMovingAverages(x,y,movingWindHalf,binCenters)
%creates moving window average approximation of x,y dataset
minNumPtsAvg=5;
minNumPtsAvg=2;
%minNumPtsAvg=1;

numApproxPts=length(binCenters);
yMovAvgValues=NaN(size(binCenters));
yMovSEMValues=NaN(size(binCenters));
for bi=1:numApproxPts
    currBinStart=binCenters(bi)-movingWindHalf;
    currBinEnd=binCenters(bi)+movingWindHalf;
    
    idxesInCurrBin=x>=currBinStart & x<currBinEnd;

    yInCurrBin=y(idxesInCurrBin);
  
    
    currBinYavg=nanmean(yInCurrBin);
    currBinYsem=getSEMacrossRows(yInCurrBin(:));
      %currBinYsem=nanstd(yInCurrBin(:));
    
    
     %currBinYavg=nanmedian(yInCurrBin);
     if(length(yInCurrBin(~isnan(yInCurrBin)))>=minNumPtsAvg)
        yMovAvgValues(bi)=currBinYavg;
        yMovSEMValues(bi)=currBinYsem;
     end
end




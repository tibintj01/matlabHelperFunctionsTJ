function [smoothedTimeSeries] = smoothHalfWind(timeSeries,halfWindIdx)
    %phase time series in degrees
    numPts=length(timeSeries);
    
    smoothedTimeSeries=NaN(numPts,1);
    for i=1:numPts
        currStartIdx=max(1,i-halfWindIdx);
        currEndIdx=min(i+halfWindIdx,numPts);
        currValues=timeSeries(currStartIdx:currEndIdx);
        
        %not empty or just NaNs
        if(~isempty(currValues) && sum(isnan(currValues))~=length(currValues))
            smoothedTimeSeries(i)=nanmean(currValues);
        end
    
    end
end


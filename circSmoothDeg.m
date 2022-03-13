function [smoothedPhaseTimeSeries] = circSmoothDeg(phaseTimeSeries,halfWindIdx)
    %phase time series in degrees
    numPts=length(phaseTimeSeries);
    
    smoothedPhaseTimeSeries=NaN(numPts,1);
    for i=1:numPts
        currStartIdx=max(1,i-halfWindIdx);
        currEndIdx=min(i+halfWindIdx,numPts);
        currPhases=phaseTimeSeries(currStartIdx:currEndIdx);
        
        %not empty or just NaNs
        if(~isempty(currPhases) && sum(isnan(currPhases))~=length(currPhases))
            smoothedPhaseTimeSeries(i)=circMeanDeg(currPhases);
        end
    
    end


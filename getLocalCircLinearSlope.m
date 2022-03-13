function [localCircLinearSlope]=getLocalCircLinearSlope(phaseTimeSeries,halfWindIdx)

%phase time series in degrees
    numPts=length(phaseTimeSeries);
    
    localCircLinearSlope=NaN(numPts,1);
    for i=1:numPts
        currStartIdx=max(1,i-halfWindIdx);
        currEndIdx=min(i+halfWindIdx,numPts);
        currPhases=phaseTimeSeries(currStartIdx:currEndIdx);
        
        currPhases=currPhases(~isnan(currPhases));
        %not empty or just NaNs
        if(length(currPhases)>=2)
            if(range(currPhases)>180)
                continue
            end
            if(var(currPhases)==0)
                s=0;
            else
                  coefficients = polyfit((1:length(currPhases))', currPhases, 1);
                  s=coefficients(1);
             disp('')
             
                  %[ rho,p,s,b ] = kempter_lincirc( (1:length(currPhases))', ang2rad(currPhases));
            end
            %if(p<0.05)
                localCircLinearSlope(i)=s*(2*halfWindIdx+1);
            %end
                
        end
    
    end



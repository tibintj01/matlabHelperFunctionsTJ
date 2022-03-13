function [timeAxis,gaussRateOverTime] = localGaussRateSpikeTrain(spikeTimes,gaussKernelWidthSec,minTime,maxTime)
%UNTITLED3 Summary of this function goes here
        
        timeWindow= maxTime-minTime;
        N=1000000;
        
        timeStep=timeWindow/N;
        timeAxis=((minTime+timeStep/2):timeStep:(maxTime-timeStep/2))';
        
        windSizeIdx=gaussKernelWidthSec/timeStep;
        %stepSizeIdx=stepSizeSec/timeStep;
        
        %localSpikeIdxes=round(N*spikeTimes);
        localSpikeTrain=rasterize(spikeTimes,1/timeStep,timeWindow);
        
        
         w=gausswin(ceil(windSizeIdx));
         w=w-nanmin(w);
         w=w/nanmax(w)
         w=w(:);
                
        %gaussRateOverTime=nanconv(localSpikeTrain,w,'edge');
        gaussRateOverTime=zeros(size(timeAxis));
        
        for i=1:length(timeAxis)
            if(localSpikeTrain(i)==1)
                startIdx=max(1,ceil(i-windSizeIdx/2));
                endIdx=min(length(timeAxis),floor(i+windSizeIdx/2));
                
                localGaussIdxes=(startIdx:endIdx)';
                gaussRateOverTime(localGaussIdxes)=gaussRateOverTime(localGaussIdxes)+w(1:(endIdx-startIdx+1));
            end
                
        end
        %figure; plot(timeAxis,gaussRateOverTime)
   
end


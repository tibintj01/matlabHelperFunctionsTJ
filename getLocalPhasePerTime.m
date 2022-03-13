function [localAvgPhasePerTime] = getLocalPhasePerTime(allSpikeTimesAndPhases,timeBinCenters,localWindowWidth)
%gets local running circular avg phase matching timebins in timeBinCenters
        numTimeBins=length(timeBinCenters);

        approxTimeStep=median(diff(timeBinCenters));
        spikeTimes=allSpikeTimesAndPhases(:,1);
        spikePhases=allSpikeTimesAndPhases(:,2);

        localWindowHalfWidthIdx=round((localWindowWidth/approxTimeStep)/2);

        localAvgPhasePerTime=NaN(size(timeBinCenters));
        for ti=1:numTimeBins
                localWindStartIdx=max(1,ti-localWindowHalfWidthIdx);
                localWindEndIdx=min(numTimeBins,ti+localWindowHalfWidthIdx);

                localWindStartTime=timeBinCenters(localWindStartIdx);
                localWindEndTime=timeBinCenters(localWindEndIdx);
                localWindWidth=localWindEndTime-localWindStartTime;
                
                localSpikePhases=spikePhases(spikeTimes>=localWindStartTime & spikeTimes <localWindEndTime);

                localSpikePhases(isnan(localSpikePhases))=[];
                localAvgPhase=circMeanDeg(localSpikePhases);
        
                
                localAvgPhasePerTime(ti)=localAvgPhase;
        end

end


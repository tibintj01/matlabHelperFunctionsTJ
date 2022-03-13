function [timeAxis,posAxis,gaussRateOverST] = local2DGaussRateSpikeTrain(spikeTimes,spikePositions,stLimits,gaussKernelWidthSec,gaussKernelWidthM)
%UNTITLED3 Summary of this function goes here
        
        minTime=stLimits.minTime;
        maxTime=stLimits.maxTime;
        minPos=stLimits.minPos;
        maxPos=stLimits.maxPos;
        
        timeLimit= maxTime-minTime;
        posLimit=maxPos-minPos;
        
        %N=1000000;
        %Ntime=2000000;
        Ntime=300000;
        %Ntime=2000000;
         %Ntime=1000000;
        %Npos=500;  
        %Npos=100; 
        Npos=50; 
         %  Npos=75;
           %Npos=60; 
        
        timeStep=timeLimit/Ntime;
        timeAxis=((minTime+timeStep/2):timeStep:(maxTime-timeStep/2))';
        
        posStep=posLimit/Npos;
        posAxis=((minPos-posStep/2):posStep:(maxPos+posStep/2))';
        
        Npos=length(posAxis)
        
        windSizeIdxTime=gaussKernelWidthSec/timeStep;
        windSizeIdxPos=gaussKernelWidthM/posStep;
        
        %stepSizeIdx=stepSizeSec/timeStep;
        
        %localSpikeIdxes=round(N*spikeTimes);
        localSpikeTrain=rasterize(spikeTimes,1/timeStep,timeLimit);
        
        spikeIdxes=find(localSpikeTrain);
        
        positionSpikeTrain=NaN(size(localSpikeTrain));
        
        for si=1:length(spikeIdxes)
            currSpikeIdx=spikeIdxes(si);
            positionSpikeTrain(currSpikeIdx)=spikePositions(si);
        end
         %localSpikeTrainPos=rasterize(spikePositions,1/posStep,spaceLimit);
        
        %{
         w=gausswin(ceil(windSizeIdx));
         w=w-nanmin(w);
         w=w/nanmax(w);
         w=w(:);
        %}
        
        %gaussRateOverTime=nanconv(localSpikeTrain,w,'edge');
        %gaussRateOverST=zeros(length(timeAxis),length(posAxis));
        %gaussRateOverST=localSpikeTrain(:)*localSpikeTrainPos(:)';
        
        bwST=zeros(Npos,Ntime);
        
        for ti=1:length(timeAxis)
            if(~isnan(positionSpikeTrain(ti)))
                approxPosIdx=interp1(posAxis,1:length(posAxis),positionSpikeTrain(ti));
                posIdx=round(approxPosIdx);
                bwST(posIdx,ti)=1;
            end
        end
        
        disp('convolving image...')
        tic
        posSigma=round(gaussKernelWidthM/posStep);
        timeSigma=round(gaussKernelWidthSec/timeStep);
        
        if(posSigma==0)
            posSigma=1;
        end
        gaussRateOverST=imgaussfilt(bwST,[posSigma timeSigma]);
        toc
        %{
        figure; omarPcolor(timeAxis,posAxis,gaussSpikesInST)
        ylabel('Position in field (m)')
        xlabel('Time in session (sec)')
        %}
        %figure; plot(timeAxis,gaussRateOverTime)
        
        %gaussRateOverST=[];
end


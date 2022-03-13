function [xcorrInfo] = getXcorrInfo(unit1Struct,unit2Struct,positionTimeAxis,fH)
%given spike times and spike positions for overlapping fields,
%saves cross correlogram information in struct

useGoodSpeedSpikesOnly=1;
minSpeed=0.05; %m
highSpeedOnly=0;
lowSpeedOnly=0;

if(~useGoodSpeedSpikesOnly)
    unit1FieldSpikeTimes=unit1Struct.inFieldSpikeTimes;
    unit1FieldSpikePositions=unit1Struct.inFieldSpikePositions;
    unit2FieldSpikeTimes=unit2Struct.inFieldSpikeTimes;
    unit2FieldSpikePositions=unit2Struct.inFieldSpikePositions;
else
    unit1FieldSpikeTimes=unit1Struct.inFieldSpikeTimesGoodSpeed;
    unit1FieldSpikePositions=unit1Struct.inFieldSpikePositionsGoodSpeed;
    unit2FieldSpikeTimes=unit2Struct.inFieldSpikeTimesGoodSpeed;
    unit2FieldSpikePositions=unit2Struct.inFieldSpikePositionsGoodSpeed;
end
    
  
    allSpikeSpeeds=[unit1Struct.inFieldSpikeSpeeds(:) ; unit2Struct.inFieldSpikeSpeeds(:)];
    
    allGoodSpikeSpeeds=allSpikeSpeeds;
    allGoodSpikeSpeeds(allGoodSpikeSpeeds<minSpeed)=NaN;
    
    %firstQuartSpeed=prctile(allSpikeSpeeds,25);
    %thirdQuartSpeed=prctile(allSpikeSpeeds,75);
    firstQuartSpeed=prctile(allGoodSpikeSpeeds,25);
    thirdQuartSpeed=prctile(allGoodSpikeSpeeds,75);
        %{
    firstQuartSpeed=prctile(allGoodSpikeSpeeds,50);
    thirdQuartSpeed=prctile(allGoodSpikeSpeeds,50);

    if(lowSpeedOnly)
        unit1FieldSpikeTimes(allGoodSpikeSpeeds>firstQuartSpeed)=NaN;
        unit2FieldSpikeTimes(allGoodSpikeSpeeds>firstQuartSpeed)=NaN;
        unit1FieldSpikePositions(allGoodSpikeSpeeds>firstQuartSpeed)=NaN;
        unit2FieldSpikePositions(allGoodSpikeSpeeds>firstQuartSpeed)=NaN;
    elseif(highSpeedOnly)
          unit1FieldSpikeTimes(allGoodSpikeSpeeds<thirdQuartSpeed)=NaN;
        unit2FieldSpikeTimes(allGoodSpikeSpeeds<thirdQuartSpeed)=NaN;
          unit1FieldSpikePositions(allGoodSpikeSpeeds<thirdQuartSpeed)=NaN;
        unit2FieldSpikePositions(allGoodSpikeSpeeds<thirdQuartSpeed)=NaN;
    end
    %}
    
    minTime=min(positionTimeAxis);
    maxTime=max(positionTimeAxis);
     minPos=min([unit1FieldSpikePositions(:); unit2FieldSpikePositions(:)]);
    maxPos=max([unit1FieldSpikePositions(:); unit2FieldSpikePositions(:)]);
    
    stLimits.minTime=minTime;
    stLimits.maxTime=maxTime;
    stLimits.minPos=minPos;
    stLimits.maxPos=maxPos;
    
    %gaussKernelWidthSec=0.1; %sec
    %gaussKernelWidthSec=0.075; %sec
     gaussKernelWidthSec=0.015; %sec
      gaussKernelWidthSec=0.01; %sec
    %gaussKernelWidthM=0.05; %m
    %gaussKernelWidthM=0.01; %m
    gaussKernelWidthM=0.015; %m
    gaussKernelWidthM=0.01; %m
    
    disp('getting gaussian convolution with spike times...')
   tic
    %[timeAxis,unit1GaussRateOverTime] = localGaussRateSpikeTrain(unit1FieldSpikeTimes-minTime,gaussKernelWidthSec,0,maxTime-minTime);
    %[timeAxis,unit2GaussRateOverTime] = localGaussRateSpikeTrain(unit2FieldSpikeTimes-minTime,gaussKernelWidthSec,0,maxTime-minTime);
    inputSpikeTimes1=unit1FieldSpikeTimes-minTime;
    inputSpikeTimes2=unit2FieldSpikeTimes-minTime;
    inputSpikePos1=unit1FieldSpikePositions;
    inputSpikePos2=unit2FieldSpikePositions;
    
    %{
    inputSpikeTimes1=inputSpikeTimes1(~isnan(inputSpikeTimes1));
    inputSpikeTimes2=inputSpikeTimes2(~isnan(inputSpikeTimes2));
    inputSpikePos1=inputSpikePos1(~isnan(inputSpikePos1));
    inputSpikePos2=inputSpikePos2(~isnan(inputSpikePos2));
    %}
    
    [timeAxis,posAxis,unit1GaussRateOverSpaceTime] = local2DGaussRateSpikeTrain(inputSpikeTimes1,inputSpikePos1,stLimits,gaussKernelWidthSec,gaussKernelWidthM);
    
    [timeAxis,posAxis,unit2GaussRateOverSpaceTime] = local2DGaussRateSpikeTrain(inputSpikeTimes2,inputSpikePos2,stLimits,gaussKernelWidthSec,gaussKernelWidthM);

        
    toc
    
    
    approxTimeStep=median(diff(timeAxis));
    approxPosStep=median(diff(posAxis));
    
    maxlag=round(0.5/approxTimeStep);
     maxlag=round(0.75/approxTimeStep);
     
     unit2GaussRateOverSpaceTimeNorm=unit2GaussRateOverSpaceTime-mean(unit2GaussRateOverSpaceTime(:));
     unit1GaussRateOverSpaceTimeNorm=unit1GaussRateOverSpaceTime-mean(unit1GaussRateOverSpaceTime(:));
     
    %unit1GaussRateOverTimeNorm=unit1GaussRateOverSpaceTime/max(unit1GaussRateOverSpaceTime);
     %unit2GaussRateOverTimeNorm=unit2GaussRateOverTime/max(unit2GaussRateOverTime);
     
     disp('getting gaussian cross correlogram...')
    tic
    %[xgram,lags]=xcorr(unit2GaussRateOverTime,unit1GaussRateOverSpaceTime,maxlag,'coeff')
        [xgram]=normxcorr2(unit1GaussRateOverSpaceTimeNorm,unit2GaussRateOverSpaceTimeNorm);
        
        %{
        maxShiftT=round(0.12/approxTimeStep); %sec
        maxShiftX=round(0.1/approxPosStep); %m
       [xgram]=tryAll2DShifts(unit1GaussRateOverSpaceTimeNorm,unit2GaussRateOverSpaceTimeNorm,maxShiftX,maxShiftT);
%}
    toc
    
    
    centerPos=(size(xgram,1)+1)/2;
    centerTime=(size(xgram,2)+1)/2;
    
    
  
    posStep=median(diff(posAxis));
    timeStep=median(diff(timeAxis));
    %posHalfWindow=ceil(0.25/posStep); %m
     %timeHalfWindow=ceil(0.5/timeStep); %sec
         %posHalfWindow=ceil(0.5/posStep); %m
     %timeHalfWindow=ceil(1/timeStep); %sec
        %posHalfWindowFwd=ceil(0.75/posStep); %m
     %timeHalfWindowFwd=ceil(1.5/timeStep); %sec
     %{
     posHalfWindowFwd=ceil(0.5/posStep); %m
     timeHalfWindowFwd=ceil(0.5/timeStep); %sec
     
     % posHalfWindowBack=ceil(0.25/posStep); %m
          posHalfWindowBack=posHalfWindowFwd;
     timeHalfWindowBack=ceil(0.5/timeStep); %sec    
    % timeHalfWindowBack=timeHalfWindowFwd;
    %}
    
    %posHalfWindowFwd=ceil(0.3/posStep); %m
     %timeHalfWindowFwd=ceil(0.3/timeStep); %sec
      posHalfWindowFwd=ceil(0.5/posStep); %m
     timeHalfWindowFwd=ceil(0.5/timeStep); %sec
          posHalfWindowBack=posHalfWindowFwd;
     timeHalfWindowBack=timeHalfWindowFwd;  
     
     centerPosIdxes=(centerPos-posHalfWindowBack):(centerPos+posHalfWindowFwd);
     centerTimeIdxes=(centerTime-timeHalfWindowBack):(centerTime+timeHalfWindowFwd);
     posLags=(centerPosIdxes-centerPos)*posStep;
     timeLags=(centerTimeIdxes-centerTime)*timeStep;
     badTimeIdxes=centerTimeIdxes>size(xgram,2) | centerTimeIdxes<1;
      badPosIdxes=centerPosIdxes>size(xgram,1) | centerPosIdxes <1;
      
     centerTimeIdxes(badTimeIdxes)=[];
     timeLags(badTimeIdxes)=[];
     centerPosIdxes(badPosIdxes)=[];
     posLags(badPosIdxes)=[];
     
     %figure; imagesc(xgram')
  omarPcolor(posLags,timeLags,xgram(centerPosIdxes,centerTimeIdxes)',fH)
     cb=colorbar
     colormap(gca,jet)
     xlim([min(posLags) max(posLags)])
     ylim([min(timeLags) max(timeLags)])
     
     xlabel('spatial lag (m)')
     ylabel('time lag (sec)')
     ylabel(cb,'cross correlation')
     hold on
     plot(xlim,[0 0],'k--','LineWidth',5)
      plot([0 0],ylim,'k--','LineWidth',5)
      plot(xlim,[0.12 0.12],'k--','LineWidth',2)
      
     plot([0 firstQuartSpeed],[0 1],'--','Color',getGrayRGB(),'LineWidth',3)
      plot([0 thirdQuartSpeed],[0 1],'--','Color',getGrayRGB(),'LineWidth',3)
      
       plot([0 -firstQuartSpeed],[0 1],'--','Color',getGrayRGB(),'LineWidth',3)
      plot([0 -thirdQuartSpeed],[0 1],'--','Color',getGrayRGB(),'LineWidth',3)
     %plot([0 -0.5],[0 0.5],'--','Color',getGrayRGB(),'LineWidth',3)
 
        xlim([min(posLags) max(posLags)])
     ylim([min(timeLags) max(timeLags)])
    %xlim([centerPos-posHalfWindow centerPos+posHalfWindow])
    %ylim([centerTime-timeHalfWindow centerTime+timeHalfWindow])
    disp('')
    
    xgramCenter=xgram(centerPosIdxes,centerTimeIdxes);
    
    xcorrInfo.posLags=posLags;
    xcorrInfo.timeLags=timeLags;
    xcorrInfo.xgramCenter=xgramCenter;
    xcorrInfo.centerPosIdxes=centerPosIdxes;
    xcorrInfo.centerTimeIdxes=centerTimeIdxes;
    xcorrInfo.gaussKernelWidthSec=gaussKernelWidthSec;
    xcorrInfo.gaussKernelWidthM=gaussKernelWidthM;
    xcorrInfo.stLimits=stLimits;
    xcorrInfo.posStep=posStep;
    xcorrInfo.timeStep=timeStep;
    xcorrInfo.posHalfWindowFwd=posHalfWindowFwd;
    xcorrInfo.timeHalfWindowFwd=timeHalfWindowFwd;
    xcorrInfo.timeHalfWindowBack=timeHalfWindowBack;
    xcorrInfo.posHalfWindowBack=posHalfWindowBack;
    
    %xcorrInfo.unit1GaussRateOverSpaceTime=unit1GaussRateOverSpaceTime;
    %xcorrInfo.unit1GaussRateOverSpaceTime=unit1GaussRateOverSpaceTime;
    
    %omarPcolor(posAxis,timeAxis
    %{
    timeLags=lags*approxTimeStep*1000; %msec
    figure; plot(timeAxis,unit1GaussRateOverSpaceTime,timeAxis,unit2GaussRateOverTime);
    figure; stem(timeLags,xgram)
    xlabel('time lag (msec)')
    ylabel('X corr')
    hold on
    ylim([0 0.3])
    plot([0 0],ylim,'k--','LineWidth',3)
    %}
end


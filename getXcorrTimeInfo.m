function [xcorrInfo] = getXcorrTimeInfo(unit1Struct,unit2Struct,positionTimeAxis,refUnitData,fH)
%given spike times and spike positions for overlapping fields,
%saves cross correlogram information in struct

useAbsoluteSpeed=1;
useGoodSpeedSpikesOnly=1;
useGoodSpeedSpikesOnly=0;
minSpeed=0.05; %m

useMetaFieldSpeedsByLap=1;
prctileSpeedLapThresh=1/3;
prctileSpeedLapThresh=1/4;
prctileSpeedLapThresh=0.15;

numSpikesField1PerSpeed=NaN(1,2);
numSpikesField2PerSpeed=NaN(1,2);

meanBehavTimeDiffBySpeed=NaN(1,2);

minNumTrials=5;
minNumTrials=0;
showPlot=1;

for spi=1:2
    if(spi==1)
        lowSpeedOnly=1;
          highSpeedOnly=0;
    else
        lowSpeedOnly=0;
          highSpeedOnly=1;
    end

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
    

    if(unit1Struct.fieldPosEnd-unit1Struct.fieldPosStart<0)
        fieldDirStr='leftward'
    else
        fieldDirStr='rightward'
    end
    
        lapStartTimes=refUnitData.lapStartTimesPerDir.(fieldDirStr);
        lapStopTimes=refUnitData.lapStopTimesPerDir.(fieldDirStr);
        
        numLaps=length(lapStopTimes);
        %{
        unit1FieldSpikeLapNums=getSpikeLapNums(unit1FieldSpikeTimes,lapStartTimes,lapStopTimes);
        unit2FieldSpikeLapNums=getSpikeLapNums(unit2FieldSpikeTimes,lapStartTimes,lapStopTimes);
        %}

        if(strcmp(fieldDirStr,'leftward'))
            allSpikeSpeeds=-([unit1Struct.inFieldSpikeSpeeds(:) ; unit2Struct.inFieldSpikeSpeeds(:)]);
            unit1SpikeSpeeds=-(unit1Struct.inFieldSpikeSpeeds(:));
            unit2SpikeSpeeds=-(unit2Struct.inFieldSpikeSpeeds(:));
        else
            allSpikeSpeeds=([unit1Struct.inFieldSpikeSpeeds(:) ; unit2Struct.inFieldSpikeSpeeds(:)]);
            unit1SpikeSpeeds=(unit1Struct.inFieldSpikeSpeeds(:));
            unit2SpikeSpeeds=(unit2Struct.inFieldSpikeSpeeds(:));
        end
        
        

        allGoodSpikeSpeeds=allSpikeSpeeds;
        allGoodSpikeSpeeds(allGoodSpikeSpeeds<minSpeed)=NaN;

        unit1backwardsSpikeIdxes=unit1SpikeSpeeds<0;
        unit2backwardsSpikeIdxes=unit2SpikeSpeeds<0;
        
        unit1FieldSpikeTimes(unit1backwardsSpikeIdxes)=NaN;
        unit2FieldSpikeTimes(unit2backwardsSpikeIdxes)=NaN;
        
        originalNumSpikes1=length(unit1FieldSpikeTimes);
        originalNumSpikes2=length(unit2FieldSpikeTimes);
        
        %sum(unit1backwardsSpikeIdxes)
        %sum(unit2backwardsSpikeIdxes)
        
        %firstQuartSpeed=prctile(allSpikeSpeeds,25);
        %thirdQuartSpeed=prctile(allSpikeSpeeds,75);
        
        if(~useAbsoluteSpeed)
            firstQuartSpeed=prctile(allGoodSpikeSpeeds,25);
            thirdQuartSpeed=prctile(allGoodSpikeSpeeds,75);
        else
            firstQuartSpeed=0.35;
            thirdQuartSpeed=0.7;
            firstQuartSpeed=0.3;
            thirdQuartSpeed=0.6;
        end
        

        %firstQuartSpeed=prctile(allGoodSpikeSpeeds,50);
        %thirdQuartSpeed=prctile(allGoodSpikeSpeeds,50);
        if(useMetaFieldSpeedsByLap)
            isDuringSomeLap1=false(size(unit1FieldSpikeTimes));
            isDuringSomeLap2=false(size(unit2FieldSpikeTimes));
            allLapMetaFieldAvgSpikeSpeed=NaN(numLaps,1);
            allLapMetaFieldSpikeTimeDiff=NaN(numLaps,1);
            for li=1:numLaps
                currLapStartTime=lapStartTimes(li);
                currLapStopTime=lapStopTimes(li);
                
                currLapUnit1SpikeIdxes=unit1FieldSpikeTimes>=currLapStartTime & unit1FieldSpikeTimes<=currLapStopTime;
                currLapUnit2SpikeIdxes=unit2FieldSpikeTimes>=currLapStartTime & unit2FieldSpikeTimes<=currLapStopTime;
                isDuringSomeLap1(currLapUnit1SpikeIdxes)=true;
                isDuringSomeLap2(currLapUnit2SpikeIdxes)=true;
                
                currTraversalUnit1SpikeTimeCenter=nanmean(unit1FieldSpikeTimes(currLapUnit1SpikeIdxes));
                currTraversalUnit2SpikeTimeCenter=nanmean(unit2FieldSpikeTimes(currLapUnit2SpikeIdxes));
                
                allLapMetaFieldAvgSpikeSpeed(li)=nanmean([unit1SpikeSpeeds(currLapUnit1SpikeIdxes);unit2SpikeSpeeds(currLapUnit2SpikeIdxes)]);
                allLapMetaFieldSpikeTimeDiff(li)=currTraversalUnit2SpikeTimeCenter-currTraversalUnit1SpikeTimeCenter;
            end
            
            %figure; plot(currLapMetaFieldAvgSpikeSpeed)
  
            firstQuartMetaFieldSpeed=prctile(allLapMetaFieldAvgSpikeSpeed,prctileSpeedLapThresh*100);
            thirdQuartMetaFieldSpeed=prctile(allLapMetaFieldAvgSpikeSpeed,(1-prctileSpeedLapThresh)*100);
            
            unit1FieldSpikeTimes(~isDuringSomeLap1)=NaN; %not regulated spikes, not during any lap
            unit2FieldSpikeTimes(~isDuringSomeLap2)=NaN; %not regulated spikes, not during any lap
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %ABSOLUTE SPEED SEPARATION
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %{
            firstQuartMetaFieldSpeed=0.3;
             firstQuartMetaFieldSpeed=0.35;
            %firstQuartMetaFieldSpeed=0.1;
            %thirdQuartMetaFieldSpeed=0.6;
            thirdQuartMetaFieldSpeed=0.7;
            %}
            
            if(lowSpeedOnly)
                nullLapIdxes=allLapMetaFieldAvgSpikeSpeed>firstQuartMetaFieldSpeed;
                numLowSpeedTraversals=sum(~nullLapIdxes);
                
            elseif(highSpeedOnly)
                nullLapIdxes=allLapMetaFieldAvgSpikeSpeed<thirdQuartMetaFieldSpeed;
                numHighSpeedTraversals=sum(~nullLapIdxes);
            end
            
           meanBehavTimeDiffBySpeed(spi)=nanmean(allLapMetaFieldSpikeTimeDiff(~nullLapIdxes));
            
            for li=1:numLaps
                currLapStartTime=lapStartTimes(li);
                currLapStopTime=lapStopTimes(li);
                currLapUnit1SpikeIdxes=unit1FieldSpikeTimes>=currLapStartTime & unit1FieldSpikeTimes<=currLapStopTime;
                currLapUnit2SpikeIdxes=unit2FieldSpikeTimes>=currLapStartTime & unit2FieldSpikeTimes<=currLapStopTime;
                
                if(nullLapIdxes(li))
                    unit1FieldSpikeTimes(currLapUnit1SpikeIdxes)=NaN;
                    unit2FieldSpikeTimes(currLapUnit2SpikeIdxes)=NaN;
                end
                
            end
           
            
        else
            if(lowSpeedOnly)
                %{
                unit1FieldSpikeTimes(allGoodSpikeSpeeds>firstQuartSpeed)=NaN;
                unit2FieldSpikeTimes(allGoodSpikeSpeeds>firstQuartSpeed)=NaN;
                unit1FieldSpikePositions(allGoodSpikeSpeeds>firstQuartSpeed)=NaN;
                unit2FieldSpikePositions(allGoodSpikeSpeeds>firstQuartSpeed)=NaN;
                %}
                unit1FieldSpikeTimes(unit1SpikeSpeeds>firstQuartSpeed)=NaN;
                unit2FieldSpikeTimes(unit2SpikeSpeeds>firstQuartSpeed)=NaN;
                unit1FieldSpikePositions(unit1SpikeSpeeds>firstQuartSpeed)=NaN;
                unit2FieldSpikePositions(unit2SpikeSpeeds>firstQuartSpeed)=NaN;
            elseif(highSpeedOnly)
                %{
                  unit1FieldSpikeTimes(allGoodSpikeSpeeds<thirdQuartSpeed)=NaN;
                unit2FieldSpikeTimes(allGoodSpikeSpeeds<thirdQuartSpeed)=NaN;
                  unit1FieldSpikePositions(allGoodSpikeSpeeds<thirdQuartSpeed)=NaN;
                unit2FieldSpikePositions(allGoodSpikeSpeeds<thirdQuartSpeed)=NaN;
                %}
                 unit1FieldSpikeTimes(unit1SpikeSpeeds<thirdQuartSpeed)=NaN;
                unit2FieldSpikeTimes(unit2SpikeSpeeds<thirdQuartSpeed)=NaN;
                  unit1FieldSpikePositions(unit1SpikeSpeeds<thirdQuartSpeed)=NaN;
                unit2FieldSpikePositions(unit2SpikeSpeeds<thirdQuartSpeed)=NaN;
            end
        end

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
          gaussKernelWidthSec=0.1; %sec
          gaussKernelWidthSec=0.5; %sec
          gaussKernelWidthSec=0.3; %sec
          gaussKernelWidthSec=0.5; %sec
          %gaussKernelWidthSec=0.1; %sec
          
          %gaussKernelWidthSec=0.1; %sec
        %gaussKernelWidthM=0.05; %m
        %gaussKernelWidthM=0.01; %m
        %{
        gaussKernelWidthM=0.015; %m
        gaussKernelWidthM=0.01; %m
        %}


        unit1FieldSpikeTimes(isnan(unit1FieldSpikeTimes))=[];
        unit2FieldSpikeTimes(isnan(unit2FieldSpikeTimes))=[];
        
        numSpikesField1PerSpeed(spi)=length(unit1FieldSpikeTimes);
        numSpikesField2PerSpeed(spi)=length(unit2FieldSpikeTimes);
        
        disp('getting gaussian convolution with spike times...')
       tic
        [timeAxis,unit1GaussRateOverTime] = localGaussRateSpikeTrain(unit1FieldSpikeTimes-minTime,gaussKernelWidthSec,0,maxTime-minTime);
        [timeAxis,unit2GaussRateOverTime] = localGaussRateSpikeTrain(unit2FieldSpikeTimes-minTime,gaussKernelWidthSec,0,maxTime-minTime);
        inputSpikeTimes1=unit1FieldSpikeTimes-minTime;
        inputSpikeTimes2=unit2FieldSpikeTimes-minTime;
        toc

        approxTimeStep=median(diff(timeAxis));


         %maxlag=round(2/approxTimeStep);
         maxLagTime=3;
          maxLagTime=5;
          maxlag=round(maxLagTime/approxTimeStep);

         disp('getting gaussian cross correlogram...')
        tic
        %[xgram,lags]=xcorr(unit2GaussRateOverTime,unit1GaussRateOverTime,maxlag,'coeff');
        if(spi==1)
            [xgramLowSpeed,lags]=xcorr(unit2GaussRateOverTime,unit1GaussRateOverTime,maxlag,'coeff');
        else
             [xgramHighSpeed,lags]=xcorr(unit2GaussRateOverTime,unit1GaussRateOverTime,maxlag,'coeff');
        end

        toc
         timeLags=lags*approxTimeStep; %msec

          
            %figure; plot(timeAxis,unit1GaussRateOverSpaceTime,timeAxis,unit2GaussRateOverTime);
            if(spi==1)
                if(numLowSpeedTraversals>=minNumTrials)
                %stem(timeLags,xgramLowSpeed)
                yyaxis left
                 pHLow=plot(timeLags,xgramLowSpeed,'LineWidth',4)
                  ylabel('X corr low speed')
                  ylim([0 Inf])
                else
                    showPlot=0;
                end
            else
                if(numLowSpeedTraversals>=minNumTrials && numHighSpeedTraversals>=minNumTrials)
                    %stem(timeLags,xgramHighSpeed)
                    yyaxis right
                    pHHigh=plot(timeLags,xgramHighSpeed,'LineWidth',4)
                    ylabel('X corr high speed')
                    ylim([0 Inf])
                else
                    showPlot=0;
                end
            end
            xlim([-maxLagTime maxLagTime])
             xlim([-2 2])
            xlabel('time lag (sec)')
         
       
        hold on
        %ylim([0 0.15])

end %speed loop

        [maxValLow,maxIdxLowSpeed]=max(xgramLowSpeed(:));
        [maxValHigh,maxIdxHighSpeed]=max(xgramHighSpeed(:));
        
        maxTimeLagLowSpeed=timeLags(maxIdxLowSpeed);
        maxTimeLagHighSpeed=timeLags(maxIdxHighSpeed);
        
        maxVal=max([xgramLowSpeed(:); xgramHighSpeed(:)]);
        
       
if(showPlot)
            if(maxValHigh<maxValLow)
                yyaxis left
            else
                yyaxis right
            end
            plot([0 0],[0 maxVal],'k--','LineWidth',3)
            if(~useMetaFieldSpeedsByLap)
                if(useAbsoluteSpeed)
                    legend([pHLow pHHigh],{sprintf('speed<%.2fm/s',firstQuartSpeed),sprintf('speed>%.2fm/s',thirdQuartSpeed)},'Location','northwest')
                else
                    legend([pHLow pHHigh],{'bottom quartile speed','top quartile speed'})
                end
            else
                %legend([pHLow pHHigh],{sprintf('bottom quartile traversals (<%.2fm/s)',firstQuartMetaFieldSpeed),sprintf('top quartile traversals (>%.2fm/s)',thirdQuartMetaFieldSpeed)},'Location','northwest')
                legend([pHLow pHHigh],{sprintf('bottom %d%% traversals (<%.2fm/s,n=%d)',round(prctileSpeedLapThresh*100),firstQuartMetaFieldSpeed,numLowSpeedTraversals),sprintf('top %d%% traversals (>%.2fm/s,n=%d)',round(prctileSpeedLapThresh*100),thirdQuartMetaFieldSpeed,numHighSpeedTraversals)},'Location','northwest')

            end
end


        xcorrInfo.timeLags=timeLags;
        xcorrInfo.gaussKernelWidthSec=gaussKernelWidthSec;
        xcorrInfo.approxTimeStep=approxTimeStep;
        xcorrInfo.lags=lags;
        xcorrInfo.xgramLowSpeed=xgramLowSpeed;
        xcorrInfo.xgramHighSpeed=xgramHighSpeed;
        xcorrInfo.stLimits=stLimits;
        xcorrInfo.numSpikesField1PerSpeed=numSpikesField1PerSpeed;
        xcorrInfo.numSpikesField2PerSpeed=numSpikesField2PerSpeed;
        xcorrInfo.firstQuartMetaFieldSpeed=firstQuartMetaFieldSpeed;
        xcorrInfo.thirdQuartMetaFieldSpeed=thirdQuartMetaFieldSpeed;
        
        xcorrInfo.maxTimeLagHighSpeed=maxTimeLagHighSpeed;
        xcorrInfo.maxTimeLagLowSpeed=maxTimeLagLowSpeed;
        xcorrInfo.crossCorrMaxTimeDiff=maxTimeLagLowSpeed-maxTimeLagHighSpeed;
        
        xcorrInfo.allLapMetaFieldSpikeTimeDiff=allLapMetaFieldSpikeTimeDiff;
        xcorrInfo.meanBehavTimeDiffBySpeed=meanBehavTimeDiffBySpeed;
        xcorrInfo.meanBehavTimeDiffDiff=meanBehavTimeDiffBySpeed(1)-meanBehavTimeDiffBySpeed(2);
         %xcorrInfo.highSpeedMeanBehavTimeDiff=highSpeedMeanBehavTimeDiff;
        
        xcorrInfo.allLapMetaFieldAvgSpikeSpeed=allLapMetaFieldAvgSpikeSpeed;
        
        %%%%%%%%%%%%%%%%%%%%%%
        %newest analysis flag
        %%%%%%%%%%%%%%%%%%%%%%
        xcorrInfo.separatedStartCenterEnd=1;
        %xcorrInfo.noSpikesOutsideLaps=1;
              
end

        


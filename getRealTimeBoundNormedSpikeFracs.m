function [renormedSpikeTimeFracs]=getRealTimeBoundNormedSpikeFracs(spikePhases,spikeTimeFracsOriginal)
%returns rule-based time boundary normalized spike times in field
%first time not signicantly different from 60 degrees
%and mean is less than 240 degrees, after first 3 cycles

showPlots=1;
showPlots=0;
%phaseThresh=60;
phaseThresh=90;
phaseThresh=100;
phaseThresh=60;


    numSpikes=length(spikeTimeFracsOriginal);
    
    originalTimeCenters=unique(spikeTimeFracsOriginal);
    
    ignoreFirstNcycles=min([3 length(originalTimeCenters)]);
    
    numTimeBinsOriginal=length(originalTimeCenters);
    
    tolLevel=0.01;

    
    meanPhasePerTimeBin=NaN(numTimeBinsOriginal,1);
    upperLimitDeg=NaN(numTimeBinsOriginal,1);
    lowerLimitDeg=NaN(numTimeBinsOriginal,1);
    differentFromThresh=NaN(numTimeBinsOriginal,1);
    
    for ti=1:numTimeBinsOriginal
        currTimeBinCenter=originalTimeCenters(ti);
        
        if(ti>=15)
            disp('')
            
        end
        
        currTimeBinSpikeIdxes=(abs(spikeTimeFracsOriginal-currTimeBinCenter)<tolLevel);
        
        currTimeBinSpikePhases=spikePhases(currTimeBinSpikeIdxes);
        
        if(~isempty(currTimeBinSpikePhases))
            currTimeBinMeanPhase=circMeanDeg(currTimeBinSpikePhases);
        else
            currTimeBinMeanPhase=NaN;
            continue
        end
        
        
        %[currTimeBinMeanPhaseRad, upperLimitRad, lowerLimitRad]= circ_mean(ang2rad(currTimeBinSpikePhases));
        [h mu ul ll]=circ_mtest(ang2rad(currTimeBinSpikePhases),ang2rad(phaseThresh));
        
        [h2 mu2 ul2 ll2]=circ_mtest(ang2rad(currTimeBinSpikePhases),ang2rad(phaseThresh/2));
        
        [h3 mu3 ul3 ll3]=circ_mtest(ang2rad(currTimeBinSpikePhases),ang2rad(phaseThresh/4));
        
        % currTimeBinMeanPhase
        %upperLimitDeg(ti)=upperLimitRad*180/pi;
        %lowerLimitDeg(ti)=lowerLimitRad*180/pi;
        
       differentFromThresh(ti)=(h & h2 & h3) | currTimeBinMeanPhase>240;
                %differentFromThresh(ti)=(h ) | currTimeBinMeanPhase>240;

        
        
        %{
        if(upperLimitDeg(ti)<0)
            upperLimitDeg(ti)=upperLimitDeg(ti)+360;
        end
        if(lowerLimitDeg(ti)<0)
            lowerLimitDeg(ti)=lowerLimitDeg(ti)+360;
        end
        %}
        
        
        
        meanPhasePerTimeBin(ti)=currTimeBinMeanPhase;
        
        
    end
    
        differentFromThresh(1:ignoreFirstNcycles)=1;
    %upperErr=rad2ang(angdiff(ang2rad(meanPhasePerTimeBin),ang2rad(upperLimitDeg)));
    %lowerErr=rad2ang(angdiff(ang2rad(meanPhasePerTimeBin),ang2rad(lowerLimitDeg)));
    
    if(showPlots)
        figure; plot(spikeTimeFracsOriginal,spikePhases,'k.','MarkerSize',10)
        hold on
        plot(originalTimeCenters,meanPhasePerTimeBin,'r-','LineWidth',3)
        yyaxis right
         plot(originalTimeCenters,differentFromThresh,'b-','LineWidth',3)
    end
    
    firstIdxCrossingThresh=1;
    
    %while(isnan(meanPhasePerTimeBin(firstIdxCrossingThresh)) || meanPhasePerTimeBin(firstIdxCrossingThresh)>phaseThresh )
    %while(isnan(lowerLimitDeg(firstIdxCrossingThresh)) || lowerLimitDeg(firstIdxCrossingThresh)>phaseThresh )

    while(differentFromThresh(firstIdxCrossingThresh))
        firstIdxCrossingThresh=firstIdxCrossingThresh+1;
        
        if(firstIdxCrossingThresh==length(differentFromThresh) || isnan(differentFromThresh(firstIdxCrossingThresh)))
            break
        end
    end
    
    %one more cycle after reaching minimum phase
    if(firstIdxCrossingThresh<length(differentFromThresh)-1)
        firstIdxCrossingThresh=firstIdxCrossingThresh+2;
    end
    
    newFracBoundary=originalTimeCenters(firstIdxCrossingThresh);
    %plot([newFracBoundary newFracBoundary],ylim,'k--','LineWidth',3)
     
    %shadedErrorBar(originalTimeCenters,meanPhasePerTimeBin,[upperErr(:)'; lowerErr(:)'])
        %shadedErrorBar(originalTimeCenters,meanPhasePerTimeBin,upperErr(:))

%plot(originalTimeCenters,upperLimitDeg,'b-','LineWidth',3)
%plot(originalTimeCenters,lowerLimitDeg,'b-','LineWidth',3)
    
    renormedSpikeTimeFracs=spikeTimeFracsOriginal/newFracBoundary;
    
    close all
    
    
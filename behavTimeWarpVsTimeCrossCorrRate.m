close all; clear all; clc
dataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles';
showPlots=0;
setTightSubplots_SpaceTime
%sessionNames={'Achilles_10252013','Achilles_11012013','Buddy_06272013', 'Cicero_09102014', 'Cicero_09172014','Gatsby_08282013'};

%Cicero_09172014 has only 1 strongly precessing cell used
sessionNames={'Achilles_10252013','Achilles_11012013','Buddy_06272013', 'Cicero_09102014','Gatsby_08282013'};

halfLineHeight=0.01; %raster line height

   
    slowLags=[];
    fastLags=[];
    topSpeeds=[];
    bottomSpeeds=[];
for si=1:length(sessionNames)
    %try
    currSesName=sessionNames{si};
    trackFilePath=sprintf('%sTrackProps.mat',currSesName);
    currTrackData=load(trackFilePath);
    
    fieldPairXcorrNames=fieldnames(currTrackData);
    
 
    for pi=1:length(fieldPairXcorrNames)
        currPairName=fieldPairXcorrNames{pi};
        
        
        if(~contains(currPairName,'WithBehavTimes'))
            continue
        end
        currPairData=currTrackData.(currPairName);
        
        if(~isfield(currPairData,'separatedStartCenterEnd'))
            continue
        end
        
        currBehavTimeDiff=currPairData.meanBehavTimeDiffDiff;
        currXcorrLagDiff=currPairData.crossCorrMaxTimeDiff;
        
        currBehavSlowLag=currPairData.meanBehavTimeDiffBySpeed(1);
        currXcorrSlowLag=currPairData.maxTimeLagLowSpeed;
        
        currBehavFastLag=currPairData.meanBehavTimeDiffBySpeed(2);
        currXcorrFastLag=currPairData.maxTimeLagHighSpeed;
        
        topPrctileSpeed=currPairData.thirdQuartMetaFieldSpeed;
        bottomPrctileSpeed=currPairData.firstQuartMetaFieldSpeed;
        
        if(currXcorrSlowLag<0 || currXcorrFastLag<0) % not in proper order
            continue
        end
        
        bottomSpeeds=[bottomSpeeds;bottomPrctileSpeed];
        topSpeeds=[topSpeeds;topPrctileSpeed];
        
   
        slowLags=[slowLags; currXcorrSlowLag];
        fastLags=[fastLags; currXcorrFastLag];
        %plot(currBehavTimeDiff,currXcorrLagDiff,'k.','MarkerSize',30)
        %pSlow=plot(currBehavSlowLag,currXcorrSlowLag,'b.','MarkerSize',30)
        %hold on
         %pFast=plot(currBehavFastLag,currXcorrFastLag,'r.','MarkerSize',30)
    end
    

end

numBins=50;
%numBins=100;
numBins=40;
    speedEdges=linspace(0,1.2,numBins+1);
    speedBins=edgesToBins(speedEdges);
    
    numTrialsFast=length(topSpeeds);
    numTrialsSlow=length(bottomSpeeds);
    
    figure; 
    NslowSpeed=histcounts(bottomSpeeds,speedEdges);
    NslowSpeed=smooth(NslowSpeed,5);
    pSlowSpeed=plot(speedBins,NslowSpeed,'b-','LineWidth',4)
    hold on 
    
      NfastSpeed=histcounts(topSpeeds,speedEdges);
      NfastSpeed=smooth(NfastSpeed,5);
    pFastSpeed=plot(speedBins,NfastSpeed,'r-','LineWidth',4)
    hold on 
    
    timeEdges=linspace(0,2,numBins+1);
    timeBins=edgesToBins(timeEdges);
    
      xlabel('Avg running speed between field pair (m/s)')
      ylabel('Field pair count')
    leg=legend([pSlowSpeed pFastSpeed],{sprintf('bottom 15%% speed traversals (n=%d)',numTrialsSlow), sprintf('top 15%% speed traversals (n=%d)',numTrialsFast)})
    setFigFontTo(24)
    maxFigHalfWidth
    box off
    legend boxoff
    axis tight
    
    figure; 
    Nslow=histcounts(slowLags,timeEdges);
    Nslow=smooth(Nslow);
    pSlow=plot(timeBins,Nslow,'b-','LineWidth',4)
    hold on 
    
    Nfast=histcounts(fastLags,timeEdges);
    Nfast=smooth(Nfast);
    pFast=plot(timeBins,Nfast,'r-','LineWidth',4)
    
    %xlabel('Behavioral time separation between field centers')
    %ylabel('CCG firing rate peak lag (sec)')
    xlabel('CCG firing rate peak lag between field pair (sec)')
    ylabel('Field pair count')
    leg=legend([pSlow pFast],{sprintf('bottom 15%% speed traversals (n=%d)',numTrialsSlow), sprintf('top 15%% speed traversals (n=%d)',numTrialsFast)})
    setFigFontTo(18)
    box off
    legend boxoff
    e
    
    timeWarps=slowLags./fastLags;
    numWarpBins=30;
    warpEdges=linspace(0,10,numWarpBins+1);
    figure; Nwarps=histcounts(timeWarps,warpEdges)
    warpBins=edgesToBins(warpEdges);
    Nwarps=smooth(Nwarps);
    [~,mostCommonWarpIdx]=max(Nwarps);
    
    pW=plot(warpBins,Nwarps,'k-','LineWidth',4)
    hold on
    mostCommonWarpFact=warpBins(mostCommonWarpIdx);
    %plot([mostCommonWarpFact mostCommonWarpFact],ylim,'k--','LineWidth',5)
    xlabel('Time warp factor (CCG peak lag ratio, low vs high speed)')
    ylabel('Field pair count')
    leg=legend(pW,sprintf('n=%d field pairs',numTrialsFast))
    maxFigHalfWidth
    setFigFontTo(24)
    box off
    legend boxoff
function [inFieldDelXTP,inFieldDelXTC,inFieldDelXTS,twoCycleTimeDiffPerTime,twoCycleSpaceDiffPerTime] = plotInDelCycleXDelCycleTDelCyclePspace(allElapsedDistsInFieldCm,meanSpikeTimeInField,allPhasesInField,phaseDiffPerTimeDuringField,...
    avgFieldTime,avgFieldWidth,allLapNumsInField,cycleNumsInField,timeStampsInField,absSpeedsDuringField,cycleStartTimesPerCycle,distInFieldPerTime,timeAxis,...
    fieldWidthPerLap,timeInFieldPerLap,allCycleDursInField,outlierLapNums)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %***elapsed times are too binary - use theta cycle duration times - DONE
    %***only consider diffs between theta cycles so that measures are
    %   not artifically correlated with each other
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %dispIdxes=~isnan(allPhasesInField) & abs(phaseDiffPerTimeDuringField)>0;
    minSpeedPrecession=5; %cm/s
    %minSpeedPrecession=0; %cm/s
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %each cell puts in an amt of precession based on whole field fraction, not just local measurements
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    normalizeTimeInFieldToAvgFieldTimeWidth=0;
    normalizeDistInFieldToFieldDistWidth=0; 
    
    %normFieldTimeEachLap=1;
     normFieldTimeEachLap=0;
     normFieldTimeByAvgDuration=1;
    normFieldDistEachLap=1;
    
    c=allCycleDursInField;
    s=absSpeedsDuringField;
    
    
    x=allElapsedDistsInFieldCm;
    meanXinSameCycle=NaN(size(x));
    
    t=timeStampsInField;
    meanTinSameCycle=NaN(size(x));
    twoCycleTimeDiffPerTime=NaN(size(x));
    twoCycleSpaceDiffPerTime=NaN(size(x));
    
    p=allPhasesInField;
    originalP=p;
    
    %stepHalfWind=3;
    
     stepHalfWind=1;
     stepHalfWind=0;
      
    stepHalfWind=3;
     %stepHalfWind=5;
  %stepHalfWind=4;
  stepHalfWind=0;
    smoothP=circSmoothDeg(p,stepHalfWind);
    
    stepHalfWind=2;
     %stepHalfWind=3;
    localDeltaPperTime=getLocalCircLinearSlope(smoothP,stepHalfWind);
    
    for ti=1:length(cycleNumsInField)
        
        currLapNum=allLapNumsInField(ti);
        
        if(isnan(currLapNum))
            continue
        end
        
        %outlier lap duration (skip unusally long laps)
        if(min(abs(currLapNum-outlierLapNums))==0)
            continue
        end
        currLapTotalTime=timeInFieldPerLap(currLapNum);
        currLapTotalDist=fieldWidthPerLap(currLapNum);
        
        currCycleNum=cycleNumsInField(ti);
        sameCycleIdxes=find(cycleNumsInField==currCycleNum);
        numIdxesInSameCycle=length(sameCycleIdxes(~isnan(sameCycleIdxes)));
        
        meanXinSameCycle(ti)=nansum(x(sameCycleIdxes))/numIdxesInSameCycle;
        
        meanTinSameCycle(ti)=nansum(t(sameCycleIdxes))/numIdxesInSameCycle; %is center to center of theta cycle representation of how much time has passed?
        
   
        currCycleStartTime=cycleStartTimesPerCycle(currCycleNum);
        
        try
            nextCycleEndTime=cycleStartTimesPerCycle(currCycleNum+2);
        catch ME
            disp(ME.message)
            nextCycleEndTime=cycleStartTimesPerCycle(end);
        end
        
        nextCycleEndPos=interp1(timeAxis,distInFieldPerTime,nextCycleEndTime);
        currCycleStartPos=interp1(timeAxis,distInFieldPerTime,currCycleStartTime);
        
        twoCycleTimeDiffPerTime(ti)=nextCycleEndTime-currCycleStartTime;
        %interpolate position corresponding to each time event
        twoCycleSpaceDiffPerTime(ti)=nextCycleEndPos-currCycleStartPos;
        
        if(normFieldTimeEachLap)
            twoCycleTimeDiffPerTime(ti)=twoCycleTimeDiffPerTime(ti)/currLapTotalTime;
        elseif(normFieldTimeByAvgDuration)
            twoCycleTimeDiffPerTime(ti)=twoCycleTimeDiffPerTime(ti)/avgFieldTime;
        end
        
        if(normFieldDistEachLap)
            twoCycleSpaceDiffPerTime(ti)=twoCycleSpaceDiffPerTime(ti)/currLapTotalDist;
        end
    end
    
    %dxPerTimeStep=[NaN; diff(meanXinSameCycle)];
     dxPerTimeStep=twoCycleSpaceDiffPerTime(:);
    %dtPerTimeStep=[NaN;diff(meanTinSameCycle)];
    dtPerTimeStep=twoCycleTimeDiffPerTime(:);
    %dpPerTimeStep=[NaN;diff(p)];
     dpPerTimeStep=[NaN;angdiffDeg(smoothP)];
    
     %dpPerTimeStep=localDeltaPperTime;
     
     originalX=dxPerTimeStep;
     originalT=dtPerTimeStep;

    %initial phase should be greater than amount of precession into next
    %cycle (otherwise not part of during precession)
    cycleChangeIdxes=[NaN;diff(cycleNumsInField)];
    smoothCycleChange=smoothHalfWind(cycleChangeIdxes,stepHalfWind);
    cycleVicinityThresh=1/(2*stepHalfWind+1); %cycle within this many points
    
    %goodChangeCycleIdxes=~isnan(allPhasesInField) & allPhasesInField>abs(dpPerTimeStep) & absSpeedsDuringField>minSpeedPrecession & cycleChangeIdxes==1; %only look at (delta phase, delta t, delta x) where cycle number changes by 1;
    goodChangeCycleIdxes=~isnan(allPhasesInField) & absSpeedsDuringField>minSpeedPrecession & cycleChangeIdxes==1 & ~isnan(dxPerTimeStep); %only look at (delta phase, delta t, delta x) where cycle number changes by 1;
        %goodChangeCycleIdxes=~isnan(allPhasesInField) & absSpeedsDuringField>minSpeedPrecession & smoothCycleChange>cycleVicinityThresh & ~isnan(dxPerTimeStep); %only look at (delta phase, delta t, delta x) where cycle number changes by 1;

    %goodChangeCycleIdxes=~isnan(allPhasesInField) & absSpeedsDuringField>minSpeedPrecession & ~isnan(dxPerTimeStep); 
    %now all time points have temporal information...

    %x=allElapsedDistsInFieldCm(dispIdxes);
    
    cycleChangeIdxes(~goodChangeCycleIdxes)=NaN;
    %c=cycleChangeIdxes;
    %c(~goodChangeCycleIdxes)=NaN;
    
    %cycleNumPerCycleChange=cycleNumsInField(goodChangeCycleIdxes);
    
    %dcPerCycleChange=[NaN; diff(cycleNumPerCycleChange)];
    
    %dxPerCycleChange=dxPerTimeStep(goodChangeCycleIdxes);
    %dtPerCycleChange=dtPerTimeStep(goodChangeCycleIdxes);
    %pPerCycleChange=dpPerTimeStep(goodChangeCycleIdxes);
    %li=allLapNumsInField(goodChangeCycleIdxes);
    
    dxPerCycleChange=dxPerTimeStep;
    dtPerCycleChange=dtPerTimeStep;
    dpPerCycleChange=dpPerTimeStep;
    
    dxPerCycleChange(~goodChangeCycleIdxes)=NaN;
    dtPerCycleChange(~goodChangeCycleIdxes)=NaN;
    dpPerCycleChange(~goodChangeCycleIdxes)=NaN;
    c(~goodChangeCycleIdxes)=NaN;
    s(~goodChangeCycleIdxes)=NaN;
    originalX(~goodChangeCycleIdxes)=NaN;
    originalT(~goodChangeCycleIdxes)=NaN;
    
    %dxPerCycleChange=dxPerTimeStep(goodChangeCycleIdxes);
    %dtPerCycleChange=dtPerTimeStep(goodChangeCycleIdxes);
    %dpPerCycleChange=dpPerTimeStep(goodChangeCycleIdxes);
    %li=allLapNumsInField(goodChangeCycleIdxes);
    
    
    %dxPerCycleChange(dcPerCycleChange~=1)=NaN;
    %dtPerCycleChange(dcPerCycleChange~=1)=NaN;
    %dpPerCycleChange(dcPerCycleChange~=1)=NaN;
    
    %dcPerCycleChange=dcPerCycleChange(dcPerCycleChange==1);
   
    
    %{
    figure; plot(goodChangeCycleIdxes)
    hold on
    plot(dtPerTimeStep,'k','LineWidth',3)
    figure; 
    %plot(c,'b')
    %hold on
    %plot(diff(x)./diff(t),'k-','LineWidth',2)
    subplot(3,1,1)
    histogram(dxPerCycleChange,20)
    subplot(3,1,2)
     histogram(dtPerCycleChange,20)
     subplot(3,1,3)
     histogram(dpPerCycleChange,20)
    %}
    useDiff=1;
    if(useDiff)
        x=dxPerCycleChange;
        t=dtPerCycleChange;
        p=dpPerCycleChange;

        
         %maxPhase=0;
        %minPhase=-180;
        %minPhase=-240;
        maxPhase=360;
        minPhase=-360;
        
        %maxPhase=180;
        %minPhase=-180;
        
        %maxPhase=150;
        %minPhase=-150;
        
          % maxPhase=30;
        %minPhase=-30;
         maxPhase=0;
        minPhase=-180;
     

        %minDist=6;
        minDist=0;
        %maxDist=15;
        %maxDist=12;
         %maxDist=25;
          maxDist=Inf;

        %minTime=0.075;
        %minTime=0;
        %maxTime=0.3;
         %minTime=0.05;
        %maxTime=0.2;
         %minTime=0.1;
        %maxTime=0.18;
        
        %minTime=0.15;
        %maxTime=0.35;
        
        %minTime=0.05;
        %maxTime=0.45;
          minTime=0;
        maxTime=Inf;
         goodIdxes=x>=0 & t>=0 & p<maxPhase; %whenever phase goes down a reasonable amount, what happened in space and time?
         %goodIdxes=x>=0 & t>=0;
         goodIdxes=goodIdxes & x<=maxDist & t<=maxTime & t>=minTime & p>=minPhase;
         %goodIdxes=goodIdxes & x<=maxDist & t<=maxTime & t>=minTime; %no phase restriction
  
    end
    
   
    
    %x=x(goodIdxes);
    %t=t(goodIdxes);
    %p=p(goodIdxes);
     x(~goodIdxes)=NaN;
    t(~goodIdxes)=NaN;
    p(~goodIdxes)=NaN;
    %c(~goodIdxes)=NaN;
    %s(~goodIdxes)=NaN;
    
    %{
     figure; plot(p(~isnan(x)),'r-')
     hold on; plot(allElapsedDistsInFieldCm(~isnan(x)),'k-')
     hold on; plot(smoothP(~isnan(x)),'k-','LineWidth',3)
     hold on; plot(dpPerTimeStep(~isnan(x)),'r-','LineWidth',3)
    plot(xlim,[0 0],'k--')
    %}
    %{
     figure; plot(allPhasesInField,'r-','LineWidth',3)
     hold on; plot(allElapsedDistsInFieldCm*360/50,'k-')
     plot(allPhasesInField,'ro')
     
     %hold on; plot(allElapsedDistsInFieldCm(~isnan(x)),'k-')
     hold on; plot(smoothP,'b--','LineWidth',2)
      plot(smoothP,'bo')
      
      %plot(localDeltaPperTime,'m-','LineWidth',2)
      %plot(localDeltaPperTime,'mo')
      plot(p,'m-','LineWidth',2)
      plot(p,'mo')
      %}
     %hold on; plot(dpPerTimeStep(~isnan(x)),'r-','LineWidth',3)
    plot(xlim,[0 0],'k--')
    disp('')
    
    
    %{
    isolatedIdxes=[];
    for i=2:(length(p)-1) 
        if(isnan(p(i-1)) && isnan(p(i+1)))
            isolatedIdxes=[isolatedIdxes(:);i];
        end
    end
    if(isnan(p(2)))
        isolatedIdxes=[isolatedIdxes(:);1];
    end
    if(isnan(p(end-1)))
        isolatedIdxes=[isolatedIdxes(:);length(p)];
    end
    %p(isolatedIdxes)=NaN;
    %}
    %c=cycleNumPerCycleChange(goodIdxes);

    %{
    
    figure; plot(p)
    hold on; plot(p,'ko')
    %}
    
    inFieldDelXTP=[x(:) t(:) p(:)];
    
    inFieldDelXTC=[x(:) t(:) c(:)];
    inFieldDelXTS=[x(:) t(:) s(:)];
    
    
    %{
        
    figure;
    
      subplot(3,1,1)
      %avgFieldTime=0.9
      %avgFieldTime=0.15

    scatter(t,p,30,x,'filled')
    
    xlim([minTime maxTime])
    ylim([minPhase maxPhase])
    colormap(jet)
    cb1=colorbar;
    %ylabel(cb1,'lap number')
    ylabel(cb1,'associated distance')
    caxis([minDist maxDist])
    
      subplot(3,1,2)
      scatter(x,p,30,t,'filled')
    %plot(x,p,'k.')
     %ylim([0 360])
     xlim([minDist maxDist])
     ylim([minPhase maxPhase])
     cb2=colorbar;
     ylabel(cb2,'associated time')
     caxis([minTime maxTime])
     %caxis([0.1 0.14])
     
     subplot(3,1,3)
      scatter(x,t,30,p,'filled')
    %plot(x,p,'k.')
     ylim([minTime maxTime])
     cb2=colorbar;
     ylabel(cb2,'associated phase')
     %caxis([0 360])
       caxis([minPhase maxPhase])
       xlim([minDist maxDist])
       
       
  
       
       %distance conditional p(phase,time) distribution
       
       numDistConditions=ceil((maxDist-minDist));
       distEdges=linspace(minDist,maxDist,numDistConditions+1);
       distBinCenters=edgesToBins(distEdges);
       distBinWidth=median(diff(distBinCenters));
       numDistBins=length(distBinCenters);       

       numTimeConditions=ceil((maxTime-minTime)*50);
       timeEdges=linspace(minTime,maxTime,numTimeConditions+1);
       timeBinCenters=edgesToBins(timeEdges);
       timeBinWidth=median(diff(timeBinCenters));
       numTimeBins=length(timeBinCenters);       
       
       
       
       delXdelTdelP=NaN(numDistBins,numTimeBins);

       figure
       for di=1:numDistBins
            currBinDistStart=distBinCenters(di)-distBinWidth/2;
           currBinDistEnd=distBinCenters(di)+distBinWidth/2;
           %{
           for ti=1:numTimeBins
           currBinTimeStart=timeBinCenters(ti)-timeBinWidth/2;
           currBinTimeEnd=timeBinCenters(ti)+timeBinWidth/2;
           currBinIdxes=x>=currBinDistStart & x<currBinDistEnd & t>=currBinTimeStart & t<currBinTimeEnd;

           currBinPhases=p(currBinIdxes);
           currBinPhases=currBinPhases(currBinPhases<0);
           currBinPhases=abs(currBinPhases);

           if(isempty(currBinPhases))
               continue
           end
           currBinMeanPhase=circMeanDeg(currBinPhases);          
            delXdelTdelP(di,ti)=currBinMeanPhase;
           %}
 	   
           currBinIdxes=x>=currBinDistStart & x<currBinDistEnd;
           currBinPhases=p(currBinIdxes);
           currBinTimes=t(currBinIdxes);
           
           subplot(ceil(sqrt(numDistConditions)),ceil(sqrt(numDistConditions)),di)
           
           plot(currBinTimes,currBinPhases,'k.')
            if(useDiff)
               xlim([minTime maxTime])
               ylim([minPhase maxPhase])
            else
               xlim([min(currBinTimes) min(currBinTimes)+0.25])
               ylim([min(currBinPhases) min(currBinPhases)+100])
            end
             title(sprintf('d=%.2f to %.2f cm',currBinDistStart,currBinDistEnd))
             

            
           end

       %}
       
       %{
       figure
	omarPcolor(distBinCenters,timeBinCenters,delXdelTdelP')
	colormap(jet)
	colorbar      
    caxis([0 120])
       %}
 
          
       
    
     %{
    subplot(3,1,3)
    
    scatter3(x,t,p,30,p,'filled')
    ylim([0 1.2])
     %}
    
    disp('')
    

function [inFieldDelXTP] = plotInDelXDelTDelPspace(allElapsedDistsInFieldCm,allElapsedTimesInFieldSec,allPhasesInField,phaseDiffPerTimeDuringField,avgFieldTime,avgFieldWidth,allLapNumsInField)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %***elapsed times are too binary - use theta cycle duration times - DONE
    %***only consider diffs between theta cycles so that measures are
    %   not artifically correlated with each other
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dispIdxes=~isnan(allPhasesInField) & abs(phaseDiffPerTimeDuringField)>0;
    x=allElapsedDistsInFieldCm(dispIdxes);
    t=allElapsedTimesInFieldSec(dispIdxes);
    p=allPhasesInField(dispIdxes);
    li=allLapNumsInField(dispIdxes);
    
    useDiff=0;
    useDiff=1;
    if(useDiff)
        x=diff(x);
        t=diff(t);
        p=diff(p);

        maxPhase=0;
        minPhase=-180;
        %minPhase=-240;
        %minPhase=-270;

        %minDist=6;
        minDist=0;
        maxDist=15;

        %minTime=0.075;
        %minTime=0;
        %maxTime=0.3;
         minTime=0.1;
        maxTime=0.2;
         goodIdxes=x>=0 & t>=0 & p<=0; %whenever phase goes down a reasonable amount, what happened in space and time?
         goodIdxes=goodIdxes & x<=maxDist & t<=maxTime & t>=minTime & p>=minPhase;
    else
          maxPhase=360;
        minPhase=0;

        minDist=0;
        maxDist=avgFieldWidth;

        minTime=0;
        maxTime=avgFieldTime;
         goodIdxes=x>=0 & t>=0;
    end
    
   
    
    x=x(goodIdxes);
    t=t(goodIdxes);
    p=p(goodIdxes);
    
    inFieldDelXTP=[x(:) t(:) p(:)];

    
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
    

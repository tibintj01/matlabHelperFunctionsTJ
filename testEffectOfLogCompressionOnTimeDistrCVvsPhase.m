clear all; clc; close all
%close all; clc
x=0:0.001:1;
%x=0:0.01:1;
%x=0:0.0005:1;


phaseNoiseSigma=60;
phaseNoiseSigma=45;
phaseNoiseSigma=20;
phaseNoiseSigma=0;
%phaseNoiseSigma=10;
phaseNoiseSigma=30;
phaseNoiseSigma=45;
phaseNoiseSigma=20;
phaseNoiseSigma=30;

takeSqrt=0;


numX=length(x);

xNoiseSigma=0.15;
xNoiseSigma=0.025;
xNoiseSigma=0.08;
xNoiseSigma=0.2;
xNoiseSigma=0.1;
xNoiseSigma=0.01;
xNoiseSigma=0.3;
xNoiseSigma=0.25;
xNoiseSigma=0.2;
%xNoiseSigma=0.15;
xNoiseSigma=0;
%xNoiseSigma=0.05;
xNoiseSigma=0.05;

useCVmeasure=1;
useCVmeasure=0;

divideWidthByMean=1;

%setTightSubplots_SpaceTime
setTightSubplots_Spacious

showPlots=1;


numDraws=1000;
numDraws=2000;


numTimeBins=500;
timeEdges=linspace(0,1,numTimeBins+1);
timeBinCenters=edgesToBins(timeEdges);
timeBinWidth=median(diff(timeBinCenters));

numPhaseBins=50;
numPhaseBins=100;
numPhaseBins=360;
%numPhaseBins=120;
phaseEdges=linspace(0,360,numPhaseBins+1);
phaseEdges=linspace(45,330,numPhaseBins+1);
phaseBinCenters=edgesToBins(phaseEdges);

phaseBinWidth=median(diff(phaseBinCenters));

phaseBinColors=jet(numPhaseBins);
maxProbTimePerPhaseBin=NaN(numTimeBins,1);
distrWidth=NaN(numTimeBins,1);

fullCirc=360;

yLinMean=fullCirc*(1-x);
%yLog=logTime(x-0.025);
yLogMean=fullCirc*logTime(1-x);

ylogWithNoise=NaN(numDraws,numX);
ylinWithNoise=NaN(numDraws,numX);

figBoth=figure
subplot(2,3,4)
plot(x,yLogMean,'r-','LineWidth',5)
title('phase vs log(time with gaussian noise)')

subplot(2,3,1)
plot(x,yLinMean,'k-','LineWidth',5)

maxNumDisp=500;
maxNumDisp=50;

for i=1:numDraws
   
        yLinWithNoise(i,:)=fullCirc*((1-(x+normrnd(0,xNoiseSigma,size(x)))));
        yLogWithNoise(i,:)=fullCirc*logTime(1-(x+normrnd(0,xNoiseSigma,size(x))));
  
    
    yLinWithNoise(i,:)=mod(yLinWithNoise(i,:)+normrnd(0,phaseNoiseSigma,size(x)),fullCirc);
    yLogWithNoise(i,:)=mod(yLogWithNoise(i,:)+normrnd(0,phaseNoiseSigma,size(x)),fullCirc);
    
    
    if(i<maxNumDisp)
        subplot(2,3,1)
        plot(x,yLinWithNoise(i,:),'k.')
        hold on
        subplot(2,3,4)
        plot(x,yLogWithNoise(i,:),'r.')
        hold on
    
    end
  
end

    subplot(2,3,4)
ylim([0 360])
xlabel('Time in field')
ylabel('Phase (deg)')

subplot(2,3,1)
ylim([0 360])
xlabel('Time in field')
ylabel('Phase (deg)')

%plot distributions as a function of time in field

xOneCol=repmat(x(:),numDraws,1);

yLogWithNoise=yLogWithNoise';
yLinWithNoise=yLinWithNoise';
yLogWithNoiseOneCol=yLogWithNoise(:);
yLinWithNoiseOneCol=yLinWithNoise(:); %stacks COLUMNS


figLin=figure(11)


[pxySmooth,pxyLin]=getJointDistr(xOneCol,yLinWithNoiseOneCol,timeEdges,phaseEdges,figLin);

figLog=figure(12)
[pxySmooth,pxyLog]=getJointDistr(xOneCol,yLogWithNoiseOneCol,timeEdges,phaseEdges,figLog);
    

for mi=1:2
    if(mi==1)
        pxy=pxyLin;
    else
        pxy=pxyLog
    end
    
    for phi=1:numPhaseBins
            %currPdist=circSmooth(pxySmooth(ti,:),3);
            
            currPhaseBinStart=phaseBinCenters(phi)-phaseBinWidth/2;
            currPhaseBinEnd=phaseBinCenters(phi)+phaseBinWidth/2;

            if(mi==1)
                inBinIdxes=yLinWithNoiseOneCol>=currPhaseBinStart & yLinWithNoiseOneCol<currPhaseBinEnd;
                inBinPhases=yLinWithNoiseOneCol(inBinIdxes);
            else
                inBinIdxes=yLogWithNoiseOneCol>=currPhaseBinStart & yLogWithNoiseOneCol<currPhaseBinEnd;
                inBinPhases=yLogWithNoiseOneCol(inBinIdxes);
            end
            
            inBinFieldTimes=(xOneCol(inBinIdxes));
            
            notNaNIdxes=~isnan(inBinFieldTimes);
            
            currPhaseBinKappa=circ_kappa(ang2rad(inBinFieldTimes(notNaNIdxes)));
            
            %[currPhaseBinTimeCV] = getCoeffOfVar(1-inBinFieldTimes);
            
         
            currPdist=smooth(pxy(:,phi),5);
            
            currPdist=currPdist/sum(currPdist);

            [~,maxID]=max(currPdist);

           %[currDistrWidth,risingIntersectTime,risingIntersectProb,descendingIntersectTime,descendingIntersectProb]...
           %    =getFullWidthHalfMaxLinearDistr(timeBinCenters,currPdist);
           
           currDistrMean=(1-nanmean(inBinFieldTimes));
            %currDistrMean=timeBinCenters(maxID);
           
           [currDistrWidth,risingIntersectTime,risingIntersectProb,descendingIntersectTime,descendingIntersectProb]...
               =getFullWidthHalfMaxCircDistr(timeBinCenters*360,currPdist);
           
           currDistrWidth=currDistrWidth/360;
           
            risingIntersectTime=risingIntersectTime/360;
           descendingIntersectTime=descendingIntersectTime/360;
           
           
        if(showPlots)
            
            figure(figBoth)
            if(mi==1)
                subplot(2,3,2)
            else
                subplot(2,3,5)
            end

            plot(timeBinCenters,currPdist,'-','Color',phaseBinColors(phi,:),'LineWidth',3)

            hold on

            plot(risingIntersectTime,risingIntersectProb,'.','Color',phaseBinColors(phi,:),'MarkerSize',30)
            plot(descendingIntersectTime,descendingIntersectProb,'.','Color',phaseBinColors(phi,:),'MarkerSize',30)

            plot([risingIntersectTime descendingIntersectTime], [risingIntersectProb descendingIntersectProb],'Color',phaseBinColors(phi,:),'LineWidth',3)

             if(mi==1)
                subplot(2,3,3)
            else
                subplot(2,3,6)
             end
             
             ylim([0 1]*phaseNoiseSigma/20)
             %ylim([0 270])
             %ylim([0 0.15]*xNoiseSigma/0.04)
             %xlim([0 360])
             %ylim([100 180])
            % ylim([90 190])
            %{
               ylim([0 0.006])
               
            %}
                hold on
                
                if(isnan(currPhaseBinKappa))
                    disp('')
                    
                end
                
           %set(gca, 'yscale', 'log')
           
           
           if(divideWidthByMean)
               currDistrWidth=currDistrWidth/currDistrMean;
           end
           if(useCVmeasure)
              distrWidth(phi)=currPhaseBinTimeCV;
            else
             distrWidth(phi)=currDistrWidth;
           end
        
           if(useCVmeasure)
                plot(phaseBinCenters(phi),(currPhaseBinTimeCV),'.','Color',phaseBinColors(phi,:),'MarkerSize',40)
           else
                 plot(phaseBinCenters(phi),(currDistrWidth),'.','Color',phaseBinColors(phi,:),'MarkerSize',40)
           end

            hold on

            
        end
            
        if(useCVmeasure)
            distrWidth(phi)=currPhaseBinTimeCV;
        else
            distrWidth(phi)=currDistrWidth;
        end
             
            maxProbTimePerPhaseBin(phi)=timeBinCenters(maxID);
    end
    

end


subplot(2,3,1)
title('phase vs time with gaussian noise')
subplot(2,3,4)
title('phase vs log(time with gaussian noise)')

subplot(2,3,2)
xlabel('Phase (deg)')
           ylabel('Probability')
           colormap(jet)
           cb=colorbar;
           ylabel(cb,'Time(frac)')
                   xlim([0 1])
                   
                   title('time distr. vs phase, linear')
           
subplot(2,3,5)
xlabel('Phase (deg)')
           ylabel('Probability')
           colormap(jet)
           cb=colorbar
           ylabel(cb,'log(time) (frac)')
           xlim([0 1])
           title('time distr. vs phase, log compressed ')
           
subplot(2,3,3)
 xlabel('Theta phase (deg)')
           ylabel('Time distribution width (frac)')
           colormap(jet)
           cb=colorbar
           ylabel(cb,'Theta phase (deg)')
           
           title('time variability vs phase, linear')
           
subplot(2,3,6)
 xlabel('Time in field (frac)')
     if(useCVmeasure)
         ylabel('Time distribution coefficient of variation')
     else
         
           ylabel('Time distribution half-max full width (deg)')
     end
           colormap(jet)
           cb=colorbar
           ylabel(cb,'log(time) (frac)')
           %title('phase variability vs log(time with gaussian noise)')
            title('time variability vs phase, log compressed')
           
           setFigFontTo(18)
           

maxFig
saveas(gcf,'logVsLinearTimeModelPhaseVariability.png')




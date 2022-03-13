clear all; clc
close all; clc

maxTime=0.8;
%maxTime=1;
maxTime=0.8;
maxTime=1;
maxTime=1;
x=(0:0.001:maxTime);
%x=0:0.01:1;
%x=0:0.0005:1;

usePoissonNoise=0;

phaseNoiseSigma=60;
phaseNoiseSigma=45;
phaseNoiseSigma=30;
phaseNoiseSigma=10;
phaseNoiseSigma=10;
phaseNoiseSigma=0;
phaseNoiseSigma=10;
phaseNoiseSigma=50;
phaseNoiseSigma=0;
phaseNoiseSigma=20;

takeSqrt=0;

if(usePoissonNoise)
    xNoiseSigma=100000;
end
numX=length(x);
xNoiseSigma=0.1;
xNoiseSigma=0.2;
xNoiseSigma=0.15;
xNoiseSigma=0.1;
xNoiseSigma=0.05;

xNoiseSigma=0.5;
xNoiseSigma=0.25;
xNoiseSigma=0.225;
xNoiseSigma=0.25;

xNoiseSigma=0.1;
%xNoiseSigma=0.05;
xNoiseSigma=0.2;
xNoiseSigma=0.05;
xNoiseSigma=0.1;
xNoiseSigma=0.2;
xNoiseSigma=0.05;
xNoiseSigma=0.1;
xNoiseSigma=0.05;
xNoiseSigma=0.1;


xNoiseSigma=0.15;
xNoiseSigma=0.15;
%xNoiseSigma=0.2;
xNoiseSigma=0.1;
%xNoiseSigma=0.05;
xNoiseSigma=0.2;
useKappaVarMeasure=0;
xNoiseSigma=0.25;
xNoiseSigma=0.1;
xNoiseSigma=0.3;
xNoiseSigma=0.1;
xNoiseSigma=0.2;
xNoiseSigma=0.25;

xNoiseSigma=0.15;
xNoiseSigma=0.25;
xNoiseSigma=0.3;
xNoiseSigma=0.25;
xNoiseSigma=0.3;
xNoiseSigma=0.1;

xNoiseSigma=0.02;
xNoiseSigma=0.05;
xNoiseSigma=0.03;

xNoiseSigma=0.01;
%xNoiseSigma=0.02;
xNoiseSigma=0.005;
xNoiseSigma=0.015;
xNoiseSigma=0.01;
xNoiseSigma=0.1;
xNoiseSigma=0.2;
xNoiseSigma=0.3;
xNoiseSigma=0.2;
xNoiseSigma=0.35;
xNoiseSigma=0.2;
xNoiseSigma=0.1;
xNoiseSigma=0.15;
xNoiseSigma=0.25;

%setTightSubplots_SpaceTime
setTightSubplots_Spacious

showPlots=1;

numDraws=100;
numDraws=500;
numDraws=1000;
numDraws=2000;
%numDraws=100;

%numDraws=1000;
%numDraws=5000;

numTimeBins=10;
numTimeBins=20;

numTimeBins=50;
numTimeBins=100;
numTimeBins=10;
numTimeBins=30;
%numTimeBins=10;
numTimeBins=50;
%timeEdges=linspace(0,1,numTimeBins+1);


timeEdges=linspace(0,maxTime,numTimeBins+1);
%timeEdges=linspace(0,1.5,numTimeBins+1);
timeBinCenters=edgesToBins(timeEdges);
timeBinWidth=median(diff(timeBinCenters));

numPhaseBins=30;
numPhaseBins=90;
%numPhaseBins=120;

phaseEdges=linspace(0,360,numPhaseBins+1);
phaseBinCenters=edgesToBins(phaseEdges);

timeBinColors=jet(numTimeBins)
maxProbPhasePerTimeBin=NaN(numTimeBins,1);
distrWidth=NaN(numTimeBins,1);



fullCirc=360;
yLinMean=fullCirc*(1-x);
%yLog=logTime(x-0.025);
%yLogMean=fullCirc*logTime(1-x);
noRectify=1;
%noRectify=0;
base=1.33333;
base=1.45;
base=1.5;
base=1.45;
%base=sqrt(2);

base=1.5;
base=2;
base=1.33333;
base=1.5;
base=3;
base=2;
base=1.33333;

yLogMean=fullCirc*getLogTime((1-x)/maxTime,base,noRectify);

yLogWithNoise=NaN(numDraws,numX);
yLinWithNoise=NaN(numDraws,numX);

figBoth=figure
%{
subplot(2,3,4)
plot(x,yLogMean,'r-','LineWidth',5)
title('phase vs log(time with gaussian noise)')

subplot(2,3,1)
plot(x,yLinMean,'k-','LineWidth',5)
%}

maxNumDisp=500;
maxNumDisp=50;

for i=1:numDraws
   
    
    %yLinWithNoise(i,:)=fullCirc*((1-(x+normrnd(0,xNoiseSigma,size(x)))));
    %yLogWithNoise(i,:)=mod(fullCirc*getLogTime((1-x)/maxTime+normrnd(0,xNoiseSigma,size(x)),base,noRectify),fullCirc);
    %noiseVec=cumsum(normrnd(0,xNoiseSigma,size(x)));
    %noiseVec=noiseVec(end:-1:1);
    
    noiseVec=normrnd(0,xNoiseSigma,size(x));
     %noiseVec=lognrnd(0,xNoiseSigma,size(x))-1;
    
    yLinWithNoise(i,:)=fullCirc*((1-(x+noiseVec)));
    yLogWithNoise(i,:)=mod(fullCirc*getLogTime((1-x)/maxTime+noiseVec,base,noRectify),fullCirc);
    
    %{
     for ni=1:numX
         yLinWithNoise(i,ni)=fullCirc*((1-(x(ni)+normrnd(0,xNoiseSigma*ni/numX,1))));
        yLogWithNoise(i,ni)=mod(fullCirc*getLogTime((1-x(ni))/maxTime+normrnd(0,xNoiseSigma*ni/numX,1),base,noRectify),fullCirc);
     end
    %}
    
   
    
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
caxis([0 1e-3])

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


%figLin=figure(11)

subplot(2,3,1)
[pxySmooth,pxyLin]=getJointDistr(xOneCol,yLinWithNoiseOneCol,timeEdges,phaseEdges,figBoth);

%figLog=figure(12)
subplot(2,3,4)
[pxySmooth,pxyLog]=getJointDistr(xOneCol,yLogWithNoiseOneCol,timeEdges,phaseEdges,figBoth);
    



for mi=1:2
    if(mi==1)
        pxy=pxyLin;
    else
        pxy=pxyLog
    end
    
    for ti=1:numTimeBins
            %currPdist=circSmooth(pxySmooth(ti,:),3);
            
            currTimeBinStart=timeBinCenters(ti)-timeBinWidth/2;
            currTimeBinEnd=timeBinCenters(ti)+timeBinWidth/2;
            
            inBinIdxes=xOneCol>=currTimeBinStart & xOneCol<currTimeBinEnd;
            
            if(mi==1)
                inBinPhases=yLinWithNoiseOneCol(inBinIdxes);
            else
                inBinPhases=yLogWithNoiseOneCol(inBinIdxes);
            end
            
            notNaNIdxes=~isnan(inBinPhases);
            currTimeBinKappa=circ_kappa(ang2rad(inBinPhases(notNaNIdxes)));
            %currTimeBinKappa=circ_r(ang2rad(inBinPhases(notNaNIdxes)));

            
            
            currPdist=circSmooth(pxy(ti,:),5);
            %currPdist=circSmooth(pxy(ti,:),7);

            currPdist=currPdist/sum(currPdist);


              [~,maxID]=max(currPdist);
         [phaseBinCentersPadded,currPdistPadded,padLength]=getCircPaddedDistr(phaseBinCenters,currPdist,10);
         
         
           [currDistrWidth,risingIntersectPhase,risingIntersectProb,descendingIntersectPhase,descendingIntersectProb]...
               =getFullWidthHalfMaxCircDistr(phaseBinCenters,currPdist);
           
        if(showPlots)
            
            figure(figBoth)
            if(mi==1)
                subplot(2,3,2)
            else
                subplot(2,3,5)
            end

            plot(phaseBinCentersPadded,currPdistPadded,'-','Color',timeBinColors(ti,:),'LineWidth',3)


            
            hold on

            plot(risingIntersectPhase,risingIntersectProb,'.','Color',timeBinColors(ti,:),'MarkerSize',30)
            plot(descendingIntersectPhase,descendingIntersectProb,'.','Color',timeBinColors(ti,:),'MarkerSize',30)

            plot([risingIntersectPhase descendingIntersectPhase], [risingIntersectProb descendingIntersectProb],'Color',timeBinColors(ti,:),'LineWidth',3)

            %circMeanPhasePerSpeedBin(ti)=circMeanDeg(


        

             if(mi==1)
                subplot(2,3,3)
             
            else
                subplot(2,3,6)
               
             
             end
             %ylim([0 270])
             %ylim([0 10])
             %ylim([100 180])
            % ylim([90 190])
                hold on
                
                if(isnan(currTimeBinKappa))
                    disp('')
                    
                end
                
           %set(gca, 'yscale', 'log')
           
           if(useKappaVarMeasure)
              distrWidth(ti)=currTimeBinKappa;
            else
             distrWidth(ti)=currDistrWidth;
           
           end
        
           if(useKappaVarMeasure)
            plot(timeBinCenters(ti),sqrt(currTimeBinKappa),'.','Color',timeBinColors(ti,:),'MarkerSize',40)
           else
               if(takeSqrt)
                    plot(timeBinCenters(ti),sqrt(currDistrWidth-min(distrWidth)+1),'.','Color',timeBinColors(ti,:),'MarkerSize',40)
               else
                   plot(timeBinCenters(ti),(currDistrWidth),'.','Color',timeBinColors(ti,:),'MarkerSize',40)

               end
           end

            hold on
            
           
            
        end
            
        if(useKappaVarMeasure)
              distrWidth(ti)=currTimeBinKappa;
        else
             distrWidth(ti)=currDistrWidth;
           
        end
             
            maxProbPhasePerTimeBin(ti)=phaseBinCenters(maxID);
    end
    
           
           
           
    
end

%{
subplot(2,3,1)
title('phase vs time with gaussian noise')
subplot(2,3,4)
title('phase vs log(time with gaussian noise)')
%}

subplot(2,3,2)
xlabel('Phase (deg)')
           ylabel('Probability')
           colormap(jet)
           cb=colorbar
           ylabel(cb,'Time(frac)')
                   xlim([-Inf Inf])
                   
                   title('phase distr. vs time with gaussian noise')
           
subplot(2,3,5)
xlabel('Phase (deg)')
           ylabel('Probability')
           colormap(jet)
           cb=colorbar
           ylabel(cb,'log(time) (frac)')
           xlim([-Inf Inf])
           title('phase distr. vs log(time with gaussian noise)')
           
subplot(2,3,3)
 xlabel('Time in field (frac)')
           ylabel('Phase distribution half-max full width (deg)')
           colormap(jet)
           cb=colorbar
           ylabel(cb,'Time in field (frac)')
           ylim([0 270])
           
           ylim([0 230])
           
           title('phase variability vs time with gaussian noise')
           
subplot(2,3,6)
 xlabel('Time in field (frac)')
           ylabel('Phase distribution half-max full width (deg)')
           colormap(jet)
           cb=colorbar
           ylabel(cb,'log(time) (frac)')
           title('phase variability vs log(time with gaussian noise)')
            ylim([0 270])
           setFigFontTo(18)
           

maxFig
saveas(gcf,'logVsLinearTimeModelPhaseVariability.png')




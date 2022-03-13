close all
maxTime=1;
nPts=5000;
t=linspace(0,maxTime,nPts);
%dSlow=linspace(0,0.5,500);
d1AtMaxTime=1;
d2AtMaxTime=d1AtMaxTime*4;

dSlow=linspace(0,d1AtMaxTime,nPts);
dFast=linspace(0,d2AtMaxTime,nPts);
setPastelColors

yMaxDisp=1;
xMaxDisp=1;

minLogTdisp=0.01;
maxLogTdisp=1;

addLegend=0;

%{
d1=0.2;
d2=0.3;

d1=0.1;
d2=0.4;

d1=0.1;
d2=0.9;

d1=0.15;
d2=0.85;
%}


fieldCentersPerField=linspace(0.05, 0.95,4);

%fieldCentersPerField=linspace(0.1, 0.9,4);

%fieldColors=jet(4);
%fieldColors=[pastelNavyBlue; pastelGreen;pastelOrange;pastelMaroonRed];

%fieldColors=[pastelNavyBlue; pastelGreen;pastelGoldenYellow;pastelMaroonRed];
fieldColors=omar4ColorsColdToHot;
fieldColors=timeWarpPaperColorsColdToHot;

showLegend=0;

figure; 
subplot(1,2,1)

useLogScaleTime=0;
useHighSpeed=0;
plotDvsT
  

%subplot(2,2,3)
useLogScaleTime=0;
useHighSpeed=1;
plotDvsT

if(addLegend)
 legend([s1 s2], {'speed',sprintf('speed x %d',round(d2AtMaxTime/d1AtMaxTime))},'Location','best')
    legend boxoff
end

    
    %plot(
    
subplot(1,2,2)

useLogScaleTime=1;
useHighSpeed=0;
plotDvsT

%subplot(2,2,4)
useLogScaleTime=1;
useHighSpeed=1;
plotDvsT






setFigFontTo(12)
maxFig
saveas(gcf,'schematicDvsT_LogScaleVsLinearScale.png')


%legend([s1 s2], {'speed',sprintf('speed x %d',round(d2AtMaxTime/d1AtMaxTime))})







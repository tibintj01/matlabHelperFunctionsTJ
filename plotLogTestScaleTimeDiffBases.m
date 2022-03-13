
close all; clear all; clc
setTightSubplots_Spacious

load('speedVsOffsetOct25.mat')


xLog=0:0.001:1;



kVector=linspace(1,2,10);
%kVector=linspace(0.75,2.2,10);
%kVector=linspace(0.5,2,10);
%kVector=linspace(1.1,2,10);

%kVector=linspace(0.5,1.5,10);
numScalings=length(kVector);

base=sqrt(2);
%base=1.7;
%base=1.3;
noRectify=1;
noLimit=1;

bases=(1.2):0.01:(1.6);

numBases=length(bases);

baseColors=getBiColorMap(numBases);

useOneMinusX=1;
figure;
for bi=1:numBases
    
    base=bases(bi);





offsets=NaN(numScalings,1);

for ki=1:numScalings
   
    k=kVector(ki);
    if(useOneMinusX)
        yLog=getLogTime(min(kVector)*(1-xLog),base,noRectify,noLimit);
        yLogScaleX=getLogTime(k*(1-xLog),base,noRectify,noLimit);

    else
        yLog=getLogTime(xLog,base,noRectify,noLimit);
        yLogScaleX=getLogTime(k*(xLog),base,noRectify,noLimit);

    end

    
   


     
     currOffest=nanmedian(yLogScaleX-yLog);
     offsets(ki)=currOffest;
     
     %plot(k,currOffest,'.','Color',scalingColors(ki,:),'MarkerSize',40)
     %plot(k,currOffest,'.','Color',scalingColors(ki,:),'MarkerSize',40)
     hold on
     
     
    box off
    %ylim([0 0.35])
end

plot(kVector,offsets,'--','LineWidth',2,'color',baseColors(bi,:))
hold on

plot(dataSpeedRatioCenters,dataOffsetCycleFracs,'.','Color',getUmichGoldColor,'MarkerSize',40)


colormap(getBiColorMap)
cb=colorbar
ylabel(cb,'logarithm base')

caxis([min(bases) max(bases)])
%yLog=scaledata(yLog,0,1);

%legend('log(x)',sprintf('log(x) + log(%.1f) - 1',k),sprintf('log(%.1fx)',k))

%xlim([0 1])

xlabel('Input scaling factor (speed ratio)')
ylabel('Offset of logarithm (frac. of theta cycle)')
%axis tight
box off




end

plot(dataSpeedRatioCenters,dataOffsetCycleFracs,'.','Color',getUmichGoldColor,'MarkerSize',80)
pL=plot(dataSpeedRatioCenters,dataOffsetCycleFracs,'-','Color',getUmichGoldColor,'LineWidth',4)
ylim([0 0.35])
xlim([1 2])

totalFieldCount=392; %Oct 25, 2021 TJ
legend(pL,sprintf('offset vs speed (n=%d fields)',totalFieldCount),'Location','northwest')
legend boxoff
setFigFontTo(18)


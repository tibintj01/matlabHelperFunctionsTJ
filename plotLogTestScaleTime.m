
close all; clear all; clc
setTightSubplots_Spacious


xLog=0:0.001:1;



kVector=linspace(1,2,5);
%kVector=linspace(0.5,2,10);
%kVector=linspace(1.1,2,10);

%kVector=linspace(0.5,1.5,10);
numScalings=length(kVector);

base=sqrt(2);
%base=1.7;
%base=1.3;
noRectify=1;
noLimit=1;

useOneMinusX=1;

%scalingColors=copper(numScalings);
scalingColors=getBiColorMap(numScalings);

figure;

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

    
     subplot(1,2,1)
    %plot(xLog,yLog,'k-','LineWidth',3)
    hold on
    %plot(xLog,yLog+getLogTime(1.5,base,noRectify,noLimit)-1,'b-','LineWidth',3)
    %plot(xLog,yLogTwiceX,'r--','LineWidth',3)
    plot(xLog,yLogScaleX,'-','Color',scalingColors(ki,:),'LineWidth',3)
    hold on
    
    box off
    ylim([0 1.3])

     subplot(1,2,2)
     
     currOffest=nanmedian(yLogScaleX-yLog);
     offsets(ki)=currOffest;
     
     %plot(k,currOffest,'.','Color',scalingColors(ki,:),'MarkerSize',40)
     %plot(k,currOffest,'.','Color',scalingColors(ki,:),'MarkerSize',40)
     hold on
     
     
    box off
    %ylim([0 0.35])
end

plot(kVector,offsets,'k-','LineWidth',3)
subplot(1,2,1)
%colormap(copper)
colormap(getBiColorMap)
cb=colorbar('north')
ylabel(cb,'x scale factor')
caxis([min(kVector) max(kVector)])
%yLog=scaledata(yLog,0,1);

%legend('log(x)',sprintf('log(x) + log(%.1f) - 1',k),sprintf('log(%.1fx)',k))

%xlim([0 1])
%ylim([0 1.5])
xlabel('x (a.u.)')
ylabel('y=log(1-x) (a.u.)')
%axis tight
box off


setFigFontTo(18)


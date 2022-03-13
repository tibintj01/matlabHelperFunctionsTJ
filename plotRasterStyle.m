function [pH]=plotRasterStyle(eventTimes,rowNum,startTick,endTick,tickColor,tickLineWidth)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

numEvents=length(eventTimes);

if(isnan(tickColor))
    tickColor=[0 0 0];
end

if(~exist('tickLineWidth','var'))
    tickLineWidth=1;
end
if(~exist('startTick','var') || isnan(startTick))
    tickLength=1;
    startTick=rowNum-tickLength/2;
    endTick=rowNum+tickLength/2;
end

    for i=1:numEvents
        if(size(tickColor,1)==1)
            %plot([eventTimes(i) eventTimes(i)],[startTick endTick],'-','Color',tickColor,'LineWidth',2)
            %plot([eventTimes(i) eventTimes(i)],[startTick endTick],'-','Color',tickColor,'LineWidth',3)
             %plot([eventTimes(i) eventTimes(i)],[startTick endTick],'-','Color',tickColor,'LineWidth',4)
             pH=plot([eventTimes(i) eventTimes(i)],[startTick endTick],'-','Color',tickColor,'LineWidth',tickLineWidth)
        else
            %plot([eventTimes(i) eventTimes(i)],[startTick endTick],'-','Color',tickColor(i,:),'LineWidth',2)
            %plot([eventTimes(i) eventTimes(i)],[startTick endTick],'-','Color',tickColor(i,:),'LineWidth',3)
            pH=plot([eventTimes(i) eventTimes(i)],[startTick endTick],'-','Color',tickColor(i,:),'LineWidth',tickLineWidth)
        end
            hold on
    end

%colormap(gca,magma)


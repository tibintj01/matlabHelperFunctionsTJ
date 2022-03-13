clear all; close all; clc
logShift=6;

%x=0:0.005:1.25; y=log(x); 

xMax=1.8;
xMax=1.65;
xMax=1.66;
%xMax=1.65;
%x=0:0.005:1.5; 
x=0:0.005:xMax; 


%x=0:0.005:10; 
%xMax=max(x);
timeRescale=xMax;
%y=(log(x/timeRescale)+logShift)/logShift; 
%y=shiftedScaledLog(x,timeRescale,logShift);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%shiftedScaledLog returns:
%y=(log((timeRescale-x)/timeRescale)+logShift)/logShift; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y=shiftedScaledLog(x,timeRescale,logShift);

orangeRGB=[255,165,0]/255;
indigoRGB=[93 111 211]/255;
purpleRGB=[0.4940, 0.1840, 0.5560];

goldRGB=getUmichGoldColor();
%cellRGBs=[1 0 0; 0 0 1; 0 1 0; goldRGB];
%cellRGBs=[ 1 0 0; orangeRGB; 1 1 0;  0 1 0; 0 0 1; indigoRGB; purpleRGB; goldRGB];
cellRGBs=[goldRGB; 1 0 0; orangeRGB; 1 1 0;  0 1 0; 0 0 1; indigoRGB; purpleRGB];
%cellRGBs=flipud(cellRGBs);

numCells=size(cellRGBs,1);

initialVal=.01;
%initialVal=.00000001;
stepSize=.08;
stepSize=.03;
%stepSize=.2;
%stepSize=.02;
%stepSize=.01;
%stepSize=.1;
%stepSize=.4;

%stepSize=10;

numScalings=60;
%numScalings=10;
%scalings=linspace(1,5,numScalings);
scalings=linspace(1,16,numScalings);

%scalings=linspace(1,50,numScalings);
%scalings=linspace(0.1,30,10);
%scalings=linspace(1,18,numScalings);

numScalings=length(scalings);

seq=NaN(numCells,numScalings);

for s=1:numScalings
    for i=2:(numCells)
        seq(i,(numScalings-s)+1)=xMax- (initialVal+(i-1)*stepSize*scalings(s));
        %seq(i,(numScalings-s)+1)=(initialVal+(i-1)*stepSize*scalings(s));
    end
end
%seq1=[initialVal 1 2 3 4]/10;
%seq2=[initialVal 1 2 3 4]/10; 
seq=seq;

figure; 
%plot(x,y,'k-','LineWidth',5)

%plot(seq(:,1),0.9*shiftedScaledLog(seq(:,1),timeRescale,logShift),'k-','LineWidth',5)
hold on

 extraScaleFact=2.5;
  extraScaleFact=0.5;
  
  extraScaleFact=0.475;
  %extraScaleFact=0.4;
  %extraScaleFact=0.35;
   extraShift=-0;
 seq=seq*extraScaleFact+extraShift;
 timeRescale=timeRescale*extraScaleFact;
 logShift=logShift*extraScaleFact;
 

 for s=1:numScalings
    plot(seq(:,s)+1,0.9*shiftedScaledLog(seq(:,s),timeRescale,logShift),'k-','LineWidth',8)
        hold on
 end
 
        
for s=1:numScalings
        if(s>1)
            %prevScaling=scatter(seq(:,s-1),(log(seq(:,s-1)/timeRescale)+logShift)/logShift,1000,cellRGBs,'filled');
             %prevScaling=scatter(seq(:,s-1),shiftedScaledLog(seq(:,s-1),timeRescale,logShift),1000,cellRGBs,'filled');
            
              % alpha(prevScaling,0.2)
        end
         %currScaling=scatter(seq(:,s),(log(seq(:,s)/timeRescale)+logShift)/logShift,1000,cellRGBs,'filled');
         %currScaling=scatter(seq(:,s),shiftedScaledLog(seq(:,s),timeRescale,logShift),1000,cellRGBs,'filled');

        currScaling=scatter(seq(:,s)+1,0.9*shiftedScaledLog(seq(:,s),timeRescale,logShift),800,cellRGBs,'filled');
        
        %plot(seq(:,s)+1,0.9*shiftedScaledLog(seq(:,s),timeRescale,logShift),'k-','LineWidth',5)
        hold on
        linePlots=[];
        for cellNum=1:numCells
            
            currXVal=seq(cellNum,s)+1;
            %currYVal=(log(seq(cellNum,s)/timeRescale)+logShift)/logShift;
            currYVal=0.9*shiftedScaledLog(seq(cellNum,s),timeRescale,logShift);
            
             p1=plot([max(xlim) currXVal],[currYVal currYVal],'LineWidth',5,'Color',cellRGBs(cellNum,:));
             p2=plot([currXVal currXVal],[min(ylim) currYVal],'LineWidth',5,'Color',cellRGBs(cellNum,:));
             linePlots=[linePlots;p1];
             linePlots=[linePlots;p2];
        end
        
        %maxFig
            maxFigMukkaalWidth
            %maxFigHalfWidth
          xlabel('Time to seq. end (sec)')
         ylabel('log[time to seq. end] (frac. of theta cycle)')
        setFigFontTo(32)
        %xlim([-Inf max(x)])
        xlim([x(2) max(x)])
        xlim([0.25 max(x)*1.03])
          xlim([0.25 max(x)*1.09])
        %xlim([x(3) max(x)])
         %xlim([-Inf 1.5])
        daspect([max(x) range(ylim) 1])
         %daspect([1.5 range(ylim) 1])
        ylim([0 1.1])
          %xticklabels({1.5,1,0.5,0})
          %xticklabels({1.8:-0.2:0})
          %xticklabels({(1.8:-0.2:0)-0.2})
          xticklabels({'1.6','1.4','1.2','1.0','0.8','0.6','0.4','0.2','0.0','-0.2'})
          xticklabels({'1.4','1.2','1.0','0.8','0.6','0.4','0.2','0.0','-0.2'})
           %xticklabels({'1.2','1.0','0.8','0.6','0.4','0.2','0.0','-0.2'})
          %xtickformat('%.1f')
          
          if(exist('prevScaling','var'))
          delete(prevScaling)
          end
          %{
        if(s<10)
         saveas(gcf,sprintf('logTimingExampleFigScaling0%d.tif',s))
        else
            saveas(gcf,sprintf('logTimingExampleFigScaling%d.tif',s))
        end
          %}
        
         if(s<10)
         %saveas(gcf,sprintf('logTimingExampleFigScaling0%d.tif',s))
         print(sprintf('logTimingExampleFigScaling0%d',s),'-depsc')
        else
            %saveas(gcf,sprintf('logTimingExampleFigScaling%d.tif',s))
            print(sprintf('logTimingExampleFigScaling%d',s),'-depsc')
         end
        
         if(s>1)
            %delete(prevScaling)
              %delete(prevScaling)
         end
         if(s<numScalings)
             delete(currScaling)
             delete(linePlots)
         end
        %hold on
end
           
         


%{
plot(xlim,[log(seq2(2)) log(seq2(2))],'b-'); 
plot(xlim,[log(seq2(3)) log(seq2(3))],'b-'); 
plot(xlim,[log(seq2(4)) log(seq2(4))],'b-'); 
plot(xlim,[log(seq2(1)) log(seq2(1))],'b-');

%}

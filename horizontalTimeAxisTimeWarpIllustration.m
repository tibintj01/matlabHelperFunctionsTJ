clear all; close all; clc

figFont=24;
logShift=6;
x=0:0.005:1.5; 
x=0:0.005:1.8; 
xMax=max(x);
timeRescale=xMax;

orangeRGB=[ 255,165,0]/255;
indigoRGB=[93 111 211]/255;
purpleRGB=[0.4940, 0.1840, 0.5560];

goldRGB=getUmichGoldColor();
%cellRGBs=[1 0 0; 0 0 1; 0 1 0; goldRGB];
cellRGBs=[ 1 0 0; orangeRGB; 1 1 0;  0 1 0; 0 0 1; indigoRGB; purpleRGB; goldRGB];
cellRGBs=flipud(cellRGBs);

numCells=size(cellRGBs,1);


initialVal=.1;
stepSize=.08;
stepSize=.03;
%stepSize=.02;
stepSize=.01;

numScalings=2;
%scalings=linspace(1,5,numScalings);
scalings=linspace(1,2.5,numScalings)
%scalings=linspace(1,18,numScalings)


numScalings=length(scalings);

seq=NaN(numCells,numScalings);


for s=1:numScalings
    for i=2:(numCells)
        seq(i,(numScalings-s)+1)=xMax- (initialVal+(i-1)*stepSize*scalings(s));
    end
end

figure;
for s=1:numScalings
        initialVal=min(seq(:,s))-0.005;
         currScaling=scatter((seq(:,s)-initialVal)*3,s*ones(size(seq(:,s))),500,cellRGBs,'filled','MarkerEdgeColor','k');
         %currScaling=scatter((seq(:,s)*3),s*ones(size(seq(:,s))),500,cellRGBs,'filled','MarkerEdgeColor','k');

         hold on
         %shiftedScaledLog(seq(:,s),timeRescale,logShift)
         
end
ylim([0.1 numScalings+0.9])
curtick = get(gca, 'yTick');
yticks(unique(round(curtick)));
%xticklabels([0.6 0.5 0.4  0.3 0.2  0.1  0]*3)
%yticklabels({1 2})
xlabel('Time to seq. end (sec)')
ylabel('Sequence number (different time-warp)')

title('Feature occurence times')

maxFigHalfWidth
setFigFontTo(figFont)

figure;
for s=1:numScalings
        initialVal=min(seq(:,1))-0.005;
        rescaleVal=max((seq(:,s)-initialVal+0.1)*3)+0.1;
         %currScaling=scatter(shiftedScaledLog((seq(:,s)-initialVal+0.1)*3,rescaleVal,logShift),s*ones(size(seq(:,s))),500,cellRGBs,'filled','MarkerEdgeColor','k');
               
         flippedOrder=flipud(seq(:,s));
         flippedOrder=[NaN;flippedOrder(:)];
         flippedOrder(end)=[];
         logScaling=scatter(shiftedScaledLog(flippedOrder,timeRescale,logShift),s*ones(size(seq(:,s))),500,cellRGBs,'filled','MarkerEdgeColor','k');

         hold on
         
         
end
ylim([0.1 numScalings+0.9])
curtick = get(gca, 'yTick');
yticks(unique(round(curtick)));

%yticklabels({1 2})
xlabel('Spike timing (a.u.)')
ylabel('Sequence number (different time-warp)')

title('Feature representation log(time)')

maxFigHalfWidth
setFigFontTo(figFont)




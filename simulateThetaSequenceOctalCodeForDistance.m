close all; clear all
delX=0.1;

x=(0:delX:100);


denom=7;

%offsetFact=10;
%offsetFact=15;
offsetFact=14.275;

precessionType='log';
%precessionType='linear';

V=VideoWriter(sprintf('thetaSeqAsOctalCodes_%sPrecession.avi',precessionType));
V.FrameRate=3;
open(V)


numCells=7;
offsets=NaN(numCells,length(x));

for i=1:numCells
    offsetVal=offsetFact*(numCells-i);
    if(strcmp(precessionType,'log'))
        offsets(i,:)=(log(x-offsetVal)+2.4)*360/denom;
    else
        offsets(i,:)=3.6*(x-offsetVal);
    end
    offsets(i,x<offsetVal)=NaN;
end

moveDx=1;
%moveDx=5;
%moveDx=10;
moving=100:-moveDx:0;

numFrames=length(moving);

frames(numFrames)=struct('cdata',[],'colormap',[]);

allOdometerReadings=NaN(length(moving),numCells);

for currPosIdx=1:length(moving)
    currPos=moving(currPosIdx);
    if(exist('pH,','var'))      
         delete(pH)
    end
    %close all
    figure
    %plot(x,fieldPhaseVsDist)
    placeColors=jet(numCells)
    for cellNum=1:numCells
        plot(x,offsets(cellNum,:),'Color',placeColors(cellNum,:),'LineWidth',4)
        hold on
        [~,currPosIdxBck]=min(abs(x-currPos));
        binAtCurrPos(cellNum)=ceil(offsets(cellNum,currPosIdxBck)/(360/7));
    end
    colormap(jet)
    cb=colorbar
    caxis([1 7])
    ylabel(cb,'place cell no.')
    
    binAtCurrPos(binAtCurrPos>7)=7;
    
    binAtCurrPos(isnan(binAtCurrPos))=0;
    binAtCurrPos(isinf(binAtCurrPos))=0;
    
    seqInt=zeros(7,1);
    for gammaBin=1:7
        gammaBoundary=360/7*(gammaBin-1);
       plot(xlim,[gammaBoundary gammaBoundary],'k--','LineWidth',4)
       text(0.5,gammaBoundary+360/7/2,sprintf('8^%d',(gammaBin-1)));
       
       cellAtThisBin=find(binAtCurrPos==gammaBin);
       if(~isempty(cellAtThisBin))
            seqInt(gammaBin)=cellAtThisBin(1);
       %seqStr=[seqStr int2str(cellAtThisBin(1))];
       end
      
    end
   
    %seqStr=int2str(fliplr(seqInt))';
    allOdometerReadings(currPosIdx,:)=binAtCurrPos;
  seqStr=fliplr(int2str(binAtCurrPos));
    
    ylim([0 360])
    ylabel('Theta phase (degrees)')
    xlabel('Position (a.u.)')
    setFigFontTo(18)
    
    pH=plot([currPos currPos],ylim,'k','LineWidth',4)
    if(currPosIdx==length(moving))
        %compressedReadings
         uniqueCodeWords=unique(allOdometerReadings,'rows')
         numUniqueCodeWords=size(uniqueCodeWords,1);
         
        title({sprintf('%d unique sequences',numUniqueCodeWords),sprintf('neural sequence \"odometer\" reading:'),seqStr})
    else
     title({sprintf('%d',currPos),sprintf('neural sequence \"odometer\" reading:'),seqStr})
    end
     %set(pH,'Visible','on')
     drawnow
     
     %frames=[frames getframe(gcf)];
     %frames(currPosIdx)=getframe(gcf);
     writeVideo(V,getframe(gcf));
end
 %movie(gcf,frames(1:end-1),1,3)
 close(V)
 

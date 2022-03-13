fracPointIDs=1:1000;
numSamples=1000;

saveDir='/Users/tibinjohn/thetaSeq/code/hocFiles/Branco_2010';

posBezierYcurve=NaN(numSamples,length(fracPointIDs));
negBezierYcurve=NaN(numSamples,length(fracPointIDs));
 

for i=1:length(fracPointIDs)
    i
    fracPoint=fracPointIDs(i)/length(fracPointIDs);
    
     %X domain is 0 to 1
    %[curveX,curveY] = drawBezier([0 0 1-fracPoint 1],[1 1-fracPoint 0 0],numSamples);
   [curveX,curveY] = drawBezier([0 fracPoint 1 1],[0 0 1-fracPoint 1],numSamples);
    posBezierYcurve(:,i)=curveY;
    if(mod(i,100)==0)
        figure
    end
    plot(curveX,curveY,'b-')
    
    %[curveX,curveY] = drawBezier([0 fracPoint 1 1],[1 1 fracPoint 0],numSamples);
     [curveX,curveY] = drawBezier([0 0 1-fracPoint  1],[0 fracPoint 1 1],numSamples);
    negBezierYcurve(:,i)=curveY;
    hold on
    
    plot(curveX,curveY,'r-')
    title(sprintf('%.2f',fracPoint))
    drawnow
end

%csvwrite(fullfile(saveDir,'bezier4ptLookupPos.dat'),posBezierYcurve)

dlmwrite(fullfile(saveDir,'bezier4ptLookupPosNEW.dat'),posBezierYcurve,'delimiter',' ')
dlmwrite(fullfile(saveDir,'bezier4ptLookupNegNEW.dat'),negBezierYcurve,'delimiter',' ')

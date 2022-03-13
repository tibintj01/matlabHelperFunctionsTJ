currXlim=xlim;
currYlim=ylim;

xMax=currXlim(2);
xMin=currXlim(1);

newXTicks=[];

hold on
for mi=-100:100
    currMarkerX=0.05*mi;
    
    if(currMarkerX>=currXlim(1) && currMarkerX<=currXlim(2))
        plot([currMarkerX currMarkerX],currYlim,'b--')
        alpha 0.1
        newXTicks=[newXTicks currMarkerX];
    end
end

set(gca, 'XTick', newXTicks)
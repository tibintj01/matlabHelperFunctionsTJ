currXlim=xlim;
currYlim=ylim;

xMax=currXlim(2);
xMin=currXlim(1);

newYTicks=[];

hold on
for mi=-100:100
    currMarkerY=0.05*mi;
    
    if(currMarkerY>=currYlim(1) && currMarkerY<=currYlim(2))
        plot(currXlim,[currMarkerY currMarkerY],'k--')
        %alpha 0.1
        newYTicks=[newYTicks currMarkerY];
    end
end

set(gca, 'YTick', newYTicks)
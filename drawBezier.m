function [curveX,curveY] = drawBezier(x,y,numSamples)
    drawPlot=0;
        %drawPlot=1;
    %x and y are length 4 vector control points
    t=linspace(0,1,numSamples);
    %curveX=(1-t).^2*x(1)+2*(1-t).*t*x(2)+t.^2*x(3);
    %curveY=(1-t).^2*y(1)+2*(1-t).*t*y(2)+t.^2*y(3);
    
    curveX=(1-t).^3.*x(1) + 3*((1-t).^2).*t*x(2)+3*(1-t).*t.^2*x(3)+t.^3*x(4);
    curveY=(1-t).^3.*y(1) + 3*((1-t).^2).*t*y(2)+3*(1-t).*t.^2*y(3)+t.^3*y(4);
    
    evenlySpacedXs=t;
    correspondingYs=interp1(curveX,curveY,evenlySpacedXs);
    
    curveX=evenlySpacedXs;
    curveY=correspondingYs;
    if(drawPlot)
        figure; plot(curveX,curveY,'o')
        hold on
        for i=1:length(curveX)
           plot([ curveX(i) curveX(i)], [0 curveY(i)],'k--')
        end
        xlim([0 1])
        ylim([0 1])
    end
    
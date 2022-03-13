function y=shiftedScaledLog(x,timeRescale,logShift)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

linearApprox=1;
linearApprox=0;
    if(linearApprox)
        ySmooth=(log((timeRescale-x)/timeRescale)+logShift)/logShift; 
        %breakFracs=[1/3 2/3];
          breakFracs=[1/2 6/7];
           breakFracs=[1/2 7/8];
            breakFracs=1-[1/2 (1/2)^2 (1/2)^3 (1/2)^4 (1/2)^5];
           %breakFracs=[1/2 19/20];
        
           threePtX=[min(x)];
           for k=1:length(breakFracs)
        threePtX=[threePtX min(x)+breakFracs(k)*range(x)];
           end
           threePtX=[threePtX max(x)];
        threePtY=interp1(x,ySmooth,threePtX);
        
        threePtY(end)=0;
        figure;
        
        plot(x,ySmooth,'k')
        hold on
        plot(threePtX,threePtY,'r','LineWidth',5)
        
        hold on
        for k=1:length(threePtX)
            plot([threePtX(k) threePtX(k)],ylim,'k--','LineWidth',5)
        end
        disp('')
        
    else
    y=(log((timeRescale-x)/timeRescale)+logShift)/logShift; 
    end

end


function [rho,p,slopeDegPerXunit, offsetDeg] = getCircCorrCoeff(x,p,fH)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

x=x(:);
p=p(:);

    [ rho,p,s,b ] = kempter_lincirc(x,p*(pi/180) );
    
    slopeDegPerXunit=s*360;
    offsetDeg=rad2ang(b);
    
    xRange=range(x);
    
    if(exist('fH','var') && ~isempty(fH))
           x0=min(x)-xRange*0.2;
            y0=mod(slopeDegPerXunit*x0+offsetDeg,360);
            xf=max(x)+xRange*0.2;
            yf=mod(slopeDegPerXunit*xf+offsetDeg,360);
            
            plot([x0 xf],[y0 yf],'r--','LineWidth',3)
            
            hold on
            plot([x0 xf],[y0 yf]+360,'r--','LineWidth',3)
            plot([x0 xf],[y0 yf]-360,'r--','LineWidth',3)
    end
    

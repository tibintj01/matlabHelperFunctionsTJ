function [m,b,R,p]=getLinearFit(x,y,x0,xf,colorRGB)

x=x(:);
y=y(:);
notNaNIdx=~isnan(x) & ~isnan(y);

x=x(notNaNIdx);
y=y(notNaNIdx);
showPlots=0;
if(exist('x0','var'))
    showPlots=1;
end

if(~exist('colorRGB','var'))
    colorRGB=[0 0 0];
end
    if(length(x) <3)
        m=NaN;
        b=NaN;
        R=NaN;
        p=NaN;
        return
    end

	X=[ones(length(x),1) x(:)];

	linFit=(X.'*X)\(X.'*y(:));

	m=linFit(2);
	b=linFit(1);

	[rMatr,pMatr]=corrcoef(x(:),y(:));

	R=rMatr(2,1);
	p=pMatr(2,1);
    
    if(showPlots)
    y0=m*x0+b;
    yf=m*xf+b;
    
    
    hold on
    %plot([x0 xf],[y0 yf],'k--','LineWidth',5)
       plot([x0 xf],[y0 yf],'-','Color',colorRGB,'LineWidth',3)
    end
    
    

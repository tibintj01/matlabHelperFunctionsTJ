%dx=0.1;
%x=0:dx:1;

dt=0.01;
t=0:dt:10;
t=0:dt:20;
dt=0.02;
t=0:dt:40;
dt=0.03;
t=0:dt:60;
dt=0.1;
dt=0.05;
t=0:dt:100;

%t=0:dt:15;
%s=1;
speeds=1:0.1:2;
speeds=1:0.01:2;
speeds=1:0.001:5;
%speeds=0.2:0.001:5;
speeds=0.2:0.0005:5;
speeds=0.1:0.0005:10;
speeds=0.05:0.00025:20;
%speeds=0.2:0.0015:5;
x=max(speeds)*fliplr(t);

P=zeros(length(x),length(t));
k1=1;
k2=1;

spatialOffsets=0:5;
spatialOffsets=0:2:10;

for oi=1:length(spatialOffsets)
    downsampleStep=10;
    spatialOffset=spatialOffsets(oi);

    numSpeedPhasePlots=9;
    speedPlotCount=0;
    modVal=round(length(speeds)/numSpeedPhasePlots);
    close all
    figure
    for si=1:length(speeds)
        s=speeds(si);
        crossSection=[];
        for ti=1:length(t)
            currT=t(ti);
            xi=round(s*ti);
            currX=currT*s;
            if(xi>0 && xi<=length(x))
             %P(xi,ti)=k1*currX+k2*s*currT;
             P(xi,ti)=k1*(currX-spatialOffset)/s+k2*(currT-spatialOffset/s);
             if(P(xi,ti)<0)
                 P(xi,ti)=0;
             end
             %P(xi,ti)= P(xi,ti)*360/30;
              if(mod(si,modVal)==0)
                crossSection=[crossSection; P(xi,ti)];
              end
            end
        end
        if(mod(si,modVal)==0)
            speedPlotCount=speedPlotCount+1;
            subplot(3,3,speedPlotCount)
            plot(x(1:length(crossSection)),crossSection)
             xlabel('Distance (a.u.)')
            ylabel('Phase (a.u.)')
            %xlim([0 max(x)])
            %ylim([0 20])
             %xlim([0 10])
            %ylim([0 5])
            xlim([0 max(x)])
            ylim([0 20])
            title(sprintf('Speed=%.2f',s))
        end
    end
    uberTitle('Iso-speed lines (phase precession)')
    saveas(gcf,sprintf('isospeedLines_spacetimephaseWithOffset%.2f.tif',spatialOffset))
    %figure
    %plot3(x,t,P)
    setFigFontTo(18)
    P=P*360/max(P(:));
     figure;omarPcolor(x(1:downsampleStep:end),t(1:downsampleStep:end),P((1:downsampleStep:end),(1:downsampleStep:end)))
     %shading interp
     colormap(jet)
     cb=colorbar
     ylabel(cb,'Theta phase (degrees)')
     xlabel('Distance (a.u.)')
     ylabel('Time (a.u.)')
     title(sprintf('Spatial offset = %.2f',spatialOffset))
        setFigFontTo(18)
        saveas(gcf,sprintf('spacetimephaseWithOffset%.2f.tif',spatialOffset))
     %%
     figure
     %[~,conH]=contour(x,t,P,20);
     [~,conH]=contour(x,t,P,30);
     conH.LineWidth=3;
      colormap(jet)
      cb=colorbar
     ylabel(cb,'Theta phase (degrees)')
      xlabel('Distance (a.u.)')
     ylabel('Time (a.u.)')
     title({'Isophase lines',sprintf('Spatial offset = %.2f',spatialOffset)})
        setFigFontTo(18)
        
        saveas(gcf,sprintf('isophaseLines_spacetimephaseWithOffset%.2f.tif',spatialOffset))
        
      %%
      figure
      downsampleStep=10;
      [X,T]=meshgrid(x(1:downsampleStep:end),t(1:downsampleStep:end));
      surf(X,T,P(1:downsampleStep:end,1:downsampleStep:end))
      colormap(jet)
       xlabel('Distance (a.u.)')
     ylabel('Time (a.u.)')
     zlabel('Theta phase (degrees)')
       setFigFontTo(18)
   
    title({'Path-integration based phase precession in spacetime',sprintf('Spatial offset %d',spatialOffset)})
      %view([-1 -.5 1])
      view([-75.0894 13.6552])
      daspect([1 1 1])
     saveas(gcf,sprintf('3dsurfaceSpacetimephaseWithOffset%d.tif',spatialOffset))

end
 
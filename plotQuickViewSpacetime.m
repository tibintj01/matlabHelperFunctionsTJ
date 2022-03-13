        if(contains(currDataFilePath,'atsby'))
             bestShiftDeg=0;
         else
             bestShiftDeg=180;
         end
        
        localSpikePhasePerTimeDisp=mod(localSpikePhasePerTime-bestShiftDeg,360);
        
         
         
        %scatter(localDistInFieldPerTime,localTimeInFieldPerTime,20,localSpikePhasePerTime,'filled')
        scatter(localNormDistInFieldPerTime,localNormTimeInFieldPerTime,20,localSpikePhasePerTimeDisp,'filled')
        hold on
        colormap(jet)
        cb=colorbar
        caxis([0 360])
        %xlabel('Distance in field (m)')
        %ylabel('Time in field (sec)')
        xlabel('Distance in field (frac of avg width)')
        ylabel('Time in field (frac of avg duration)')
        ylabel(cb,'Spike theta phase (deg)')
        title(sprintf('%s, %s',removeUnderscores(fileBaseName), currFieldDirection))
        setFigFontTo(18)
        
        plot(xlim,[1 1],'k--','LineWidth',3)
         plot([1 1],ylim, 'k--','LineWidth',3)
         xlim([0 1.2])
         ylim([0 1.2])
         daspect([1 1 1])
         drawnow
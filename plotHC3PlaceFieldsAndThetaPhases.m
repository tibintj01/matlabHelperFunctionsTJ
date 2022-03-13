function [saveImgPath] = plotHC3PlaceFieldsAndThetaPhases(currUnitInfo,betterDirectionStr,saveDir)

    if(~exist('saveDir'))
          saveImage=0;
    else
      saveImage=1;
    end
    saveImgPath=NaN;
    
    positionBins=currUnitInfo.positionBins;
    firingRatePerPositionRight=currUnitInfo.firingRatePerPositionRight;
    rightSpikePositions=currUnitInfo.rightSpikePositions;
    rightSpikePhases=currUnitInfo.rightSpikePhases;
    firingRatePerPositionLeft=currUnitInfo.firingRatePerPositionLeft;
    leftSpikePositions=currUnitInfo.leftSpikePositions;
    leftSpikePhases=currUnitInfo.leftSpikePhases;

   

    
    
    speedPerTimeStep=currUnitInfo.speedPerTimeStep;
    rightwardTimeIdxes=speedPerTimeStep>0.05;
    leftwardTimeIdxes=speedPerTimeStep<-0.05;
    
    rightwardSpeedPerTimeStep=speedPerTimeStep;
    leftwardSpeedPerTimeStep=abs(speedPerTimeStep);
    
    rightwardSpeedPerTimeStep(~rightwardTimeIdxes)=NaN;
    leftwardSpeedPerTimeStep(~leftwardTimeIdxes)=NaN;
    
    
    positionPerTimeStep=currUnitInfo.positionPerTimeStep;
    
    baseName=sprintf('%s_6-12Theta_Unit%dInfo_%s',currUnitInfo.sessionName,currUnitInfo.unitIDnum,currUnitInfo.unitCellType);
    

        figH=figure;
        
        
        if(strcmp(betterDirectionStr,'betterRightward'))
            saveDirStr='rightward';
        subplot(2,1,1)
        maxFiringRate=max(firingRatePerPositionLeft);
        
        plot(positionBins,firingRatePerPositionRight,'k-','LineWidth',4)
    
        
               title('rightward lap')
               xlim([min(positionBins) max(positionBins)])
                xlabel('Position in track (m)')
                ylabel('Firing rate (Hz)')
                 box on
                 
                   yyaxis right
        plot(positionPerTimeStep,rightwardSpeedPerTimeStep,'r-','LineWidth',0.5)
        ax = gca;
        ax.YColor = 'r';
        
        ylabel('Speed (m/s)')
          %ylim([0 1.4])
          ylim([0 1])
	set(gca,'XMinorTick','on')
                 
                     subplot(2,1,2)
        plot(rightSpikePositions,rightSpikePhases,'k.')
        %scatter(rightSpikePositions,rightSpikePhases,rightCycleAmps,'k')
        
        hold on
        xlim([min(positionBins) max(positionBins)])
        
           ylim([0 360])
           
          plotLinesEvery5cm
          
            ylabel('Theta phase (deg)')
             xlabel('Position in track (m)')
             box on
             
             
    

           %box off
           
        elseif(strcmp(betterDirectionStr,'betterLeftward')) 
            saveDirStr='leftward';
           
            maxFiringRate=max(firingRatePerPositionLeft);
            
                subplot(2,1,1)
            hold on
              plot(positionBins,firingRatePerPositionLeft,'k-','LineWidth',4)
              
       xlim([min(positionBins) max(positionBins)])
        xlabel('Position in track (m)')
       ylabel('Firing rate (Hz)')
              title('leftward lap')
	set(gca,'XMinorTick','on')
              box on

              yyaxis right
        plot(positionPerTimeStep,leftwardSpeedPerTimeStep,'r-','LineWidth',0.5)
        ylabel('Speed (m/s)')
        %ylim([0 1.4])
         ylim([0 1])
           
        ax = gca;
        ax.YColor = 'r';
        
        subplot(2,1,2)
        hold on
        plot(leftSpikePositions,leftSpikePhases,'k.')
        %scatter(leftSpikePositions,leftSpikePhases,leftCycleAmps,'k')
        xlim([min(positionBins) max(positionBins)])
        ylim([0 360])
        ylabel('Theta phase (deg)')
        plotLinesEvery5cm
        xlabel('Position in track (m)')
         %box off
         box on

        end

	    set(gca,'XMinorTick','on')
    

        uberTitle(removeUnderscores(baseName))
        maxFigManual2d(3,0.8)
        setFigFontTo(12)
        
        if(saveImage)
            
           if(maxFiringRate>1)
            saveImgPath=fullfile(saveDir,sprintf('%s_PlaceFieldAndPhaseVsPos_%s.png',baseName,saveDirStr));
            %print('-r150',fullfile(saveDir,sprintf('%s_PlaceFieldAndPhaseVsPos_%s.png',baseName,saveDirStr)),'-dpng') 
            print('-r150',fullfile(saveDir,sprintf('%s_PlaceFieldAndPhaseVsPos_%s.png',baseName,saveDirStr)),'-dpng') 

           else
               saveImgPath='UNDER_ACTIVE_CELL';
           end
       
       
        else
            maxFig
        end
         %saveas(gcf,fullfile(saveDir,sprintf('%s_PlaceFieldAndPhaseVsPos_%s.tif',baseName,saveDirStr)))
           %disp('')
   %        close all
           
           %imshow(fullfile(saveDir,sprintf('%s_PlaceFieldAndPhaseVsPos_%s.png',baseName,saveDirStr)))
           %disp('')
           

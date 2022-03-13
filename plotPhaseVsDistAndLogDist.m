   if(showPlots)
            subplot(2,2,2)
            %plot(allFieldToFieldEndLogTimes,allFieldSpikePhases,'k.')

            if(plotConditionalProb)
                getJointDistrGivenX(allFieldToFieldEndLogDists,allFieldSpikePhases,logTimeFieldEdges,phaseEdges,fDistH);
            else
                getJointDistr(allFieldToFieldEndLogDists,allFieldSpikePhases,logTimeFieldEdges,phaseEdges,fDistH);
            end

            xlabel('log(dist to end of field) (normalized)')
             %xlabel('log(1-time in field)')
            ylabel('Theta phase (deg)')

            %xlim([0 360])
            %xlim([-20 380])

            %title(sprintf('Theta phase  vs log(time to end), field %d (n=%d spikes)',fi,length(currFieldSpikePhases)))


            %daspect([1 360 1])
            %caxis([0 prctile(pxySmooth(:),97.5)])
            %maxFig
            %subplot(1,2,2)


            hold on
            
   end
            x0=0;
            y0=bLogDist;
            xf=1;
            yf=mLogDist+bLogDist;

            yMid=(y0+yf)/2;



            %dataMid=circMedianDeg(allFieldSpikePhases);


            bestShift=NaN;
            shifts=-5:5;
            yMids=NaN(size(shifts));
            y0s=NaN(size(shifts));
            yFs=NaN(size(shifts));
            shiftErrors=NaN(size(shifts));
            for si=1:length(shifts)

                y0Temp=y0+360*shifts(si);
               yfTemp=yf+360*shifts(si);

               %for i=1:length(currFieldSpikePhases)
               %%%%%%%%%%%%%%%%%%%%%%%%
               %PREDICTED VALUES
               %%%%%%%%%%%%%%%%%%%%%%%%
               predictedVals=interp1([0 1],[y0Temp yfTemp],allFieldToFieldEndLogDists);
               actualVals=allFieldSpikePhases;

               errors=abs(actualVals-predictedVals);


               meanError=nanmean(errors);

               shiftErrors(si)=meanError;
               %end


               %{
                yMids(si)=(y0Temp+yfTemp)/2;
                yFs(si)=yfTemp;
                y0s(si)=y0Temp;
                %}
            end
            
           %%%%%%%%%%%%%%%%%%%%%%%%
           %PREDICTED VALUES
           %%%%%%%%%%%%%%%%%%%%%%%%
            predictedVals=interp1([0 1],[y0 yf],allFieldToFieldEndLogDists);
            predictedVals=mod(predictedVals,360);

             logDistErrors=rad2ang(angdiff(ang2rad(actualVals),ang2rad(predictedVals)));
           
             %meanErrorLogTime=nanmean(abs(errors));
             %meanErrorLogDist=nanmean((errors));
             meanErrorLogDist=circMeanDegNoMod360(logDistErrors);

            %figure; histogram(errors,30)

            %[~,bestShift]=min(abs(dataMid-yMids) + abs(360-yFs)+ abs(0-y0s));
           [~,bestShift]=min(abs(shiftErrors));


            y0=y0+360*shifts(bestShift);
               yf=yf+360*shifts(bestShift);

   if(showPlots)
            plot([x0 xf], [y0 yf],'k--','LineWidth',5)
            plot([x0 xf], [y0 yf]+360,'k--','LineWidth',5)
                 plot([x0 xf], [y0 yf]-360,'k--','LineWidth',5)
                 
            title({sprintf('Theta phase vs log(dist to end), all fields (n=%d spikes)',length(allFieldSpikePhases)),sprintf('circR=%.2f, slope=%.2f deg/field, mean circ error = %.2f deg',RLogDist,mLogDist,meanErrorLogDist)})

            
             %xlabel('log(1- (time in field))')
             xlabel('log(dist in field)')
            ylabel('Theta phase (deg)')

            ylim([0 360])
            xlim([timeMinLog 1])
             caxis(cProbLims)
            axis square
            %daspect([0 360 1])
             box off
             
             subplot(2,2,1)
            %plot(allFieldToFieldEndTimes,allFieldSpikePhases,'k.')
            if(plotConditionalProb)
                getJointDistrGivenX(allFieldToFieldEndDists,allFieldSpikePhases,timeFieldEdges,phaseEdges,fDistH);
            else
                getJointDistr(allFieldToFieldEndDists,allFieldSpikePhases,timeFieldEdges,phaseEdges,fDistH);
            end

             xlabel('Dist in field (normalized)')
            ylabel('Theta phase (deg)')
   end

             x0=0;
            y0=bBehavDist;
            xf=1;
            yf=mBehavDist+bBehavDist;

             bestShift=NaN;
            shifts=-5:5;
            yMids=NaN(size(shifts));
            y0s=NaN(size(shifts));
            yFs=NaN(size(shifts));
            shiftErrors=NaN(size(shifts));
            
            for si=1:length(shifts)
                y0Temp=y0+360*shifts(si);
               yfTemp=yf+360*shifts(si);

               %for i=1:length(currFieldSpikePhases)
               %%%%%%%%%%%%%%%%%%%%%%%%
               %PREDICTED VALUES
               %%%%%%%%%%%%%%%%%%%%%%%%
               predictedVals=interp1([0 1],[y0Temp yfTemp],allFieldToFieldEndDists);
               actualVals=allFieldSpikePhases;

               errors=abs(actualVals-predictedVals);
               %errors=rad2ang(angdiff(ang2rad(actualVals),ang2rad(predictedVals)));

               meanError=nanmean(errors);

               shiftErrors(si)=meanError;
               %end

            end

            %[~,bestShift]=min(abs(dataMid-yMids) + abs(360-yFs)+ abs(0-y0s));
           [~,bestShift]=min(abs(shiftErrors));
           
           %%%%%%%%%%%%%%%%%%%%%%%%
           %PREDICTED VALUES
           %%%%%%%%%%%%%%%%%%%%%%%%
            predictedVals=interp1([0 1],[y0 yf],allFieldToFieldEndDists);
            predictedVals=mod(predictedVals,360);
            behavDistErrors=rad2ang(angdiff(ang2rad(actualVals),ang2rad(predictedVals)));

            %meanErrorBehavTime=nanmean(abs(errors));
            %meanErrorBehavDist=nanmean((errors));
            meanErrorBehavDist=circMeanDegNoMod360(behavDistErrors);
               %meanErrorBehavTime=circMeanDeg(errors);



            y0=y0+360*shifts(bestShift);
               yf=yf+360*shifts(bestShift);
   if(showPlots)
            hold on

             plot([x0 xf], [y0 yf],'k--','LineWidth',5)
                 plot([x0 xf], [y0 yf]+360,'k--','LineWidth',5)
                 plot([x0 xf], [y0 yf]-360,'k--','LineWidth',5)

             ylim([0 360])
            xlim([timeMin 1])
             caxis(cProbLims)
             box off
             axis square

            title({sprintf('Theta phase vs Dist in field, all fields (n=%d spikes)',length(allFieldSpikePhases)),sprintf('circR=%.2f, slope=%.2f deg/field, mean circ error = %.2f deg',RBehavDist,mBehavDist,meanErrorBehavDist)})

            
            
            subplot(2,2,3)
            
            if(plotConditionalProb)
                getJointDistrGivenX(allFieldToFieldEndDists,behavDistErrors,timeFieldEdges,phaseErrEdges,fDistH);
            else
                getJointDistr(allFieldToFieldEndDists,behavDistErrors,timeFieldEdges,phaseErrEdges,fDistH);
            end

             xlabel('dist in field')
            ylabel('Linear fit error (deg)')
             caxis(cProbLims)
              box off
             axis square
                 hline(0,'k--',3)

            title({sprintf('Linear fit error vs Dist in field, all fields (n=%d spikes)',length(allFieldSpikePhases))})

              subplot(2,2,4)

            if(plotConditionalProb)
                getJointDistrGivenX(allFieldToFieldEndLogDists,logDistErrors,timeFieldEdges,phaseErrEdges,fDistH);
            else
                getJointDistr(allFieldToFieldEndLogDists,logDistErrors,timeFieldEdges,phaseErrEdges,fDistH);
            end

            caxis(cProbLims)
    
                hline(0,'k--',3)
             xlabel('log(dist in field)')
            ylabel('Linear fit error (deg)')
              box off
             axis square

            title({sprintf('Linear fit vs log(dist in field), all fields (n=%d spikes)',length(allFieldSpikePhases))})
   end

            
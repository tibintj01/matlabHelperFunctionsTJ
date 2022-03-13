           if(showPlots)
               
                subplot(2,2,2)
                %plot(allFieldToFieldEndLogTimes,allFieldSpikePhases,'k.')

                if(plotConditionalProb)
                    getJointDistrGivenX(allFieldToFieldEndLogTimes,allFieldSpikePhases,logTimeFieldEdges,phaseEdges,fTimeH);
                else
                    getJointDistr(allFieldToFieldEndLogTimes,allFieldSpikePhases,logTimeFieldEdges,phaseEdges,fTimeH);
                end

                xlabel('log(time to end of field) (normalized)')
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
            y0=bLogTime;
            xf=1;
            yf=mLogTime+bLogTime;

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
               predictedVals=interp1([0 1],[y0Temp yfTemp],allFieldToFieldEndLogTimes);
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
            predictedVals=interp1([0 1],[y0 yf],allFieldToFieldEndLogTimes);
            predictedVals=mod(predictedVals,360);
            
             logTimeErrors=rad2ang(angdiff(ang2rad(actualVals),ang2rad(predictedVals)));
           
             %meanErrorLogTime=nanmean(abs(errors));
             meanErrorLogTime=nanmean((logTimeErrors));
             %meanErrorLogTime=circMeanDeg(errors);

            %figure; histogram(errors,30)

            %[~,bestShift]=min(abs(dataMid-yMids) + abs(360-yFs)+ abs(0-y0s));
           [~,bestShift]=min(abs(shiftErrors));


            y0=y0+360*shifts(bestShift);
               yf=yf+360*shifts(bestShift);


               if(showPlots)
                    plot([x0 xf], [y0 yf],'k--','LineWidth',5)
                    plot([x0 xf], [y0 yf]+360,'k--','LineWidth',5)
                         plot([x0 xf], [y0 yf]-360,'k--','LineWidth',5)

                    title({sprintf('Theta phase vs log(time to end), all fields (n=%d spikes)',length(allFieldSpikePhases)),sprintf('circR=%.2f, slope=%.2f deg/field, mean circ error = %.2f deg',RLogTime,mLogTime,meanErrorLogTime)})


                     %xlabel('log(1- (time in field))')
                     xlabel('log(time in field)')
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
                        getJointDistrGivenX(allFieldToFieldEndTimes,allFieldSpikePhases,timeFieldEdges,phaseEdges,fTimeH);
                    else
                        getJointDistr(allFieldToFieldEndTimes,allFieldSpikePhases,timeFieldEdges,phaseEdges,fTimeH);
                    end

                     xlabel('Time in field (normalized)')
                    ylabel('Theta phase (deg)')

               end
               
             
             x0=0;
            y0=bBehavTime;
            xf=1;
            yf=mBehavTime+bBehavTime;

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
               predictedVals=interp1([0 1],[y0Temp yfTemp],allFieldToFieldEndTimes);
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
            predictedVals=interp1([0 1],[y0 yf],allFieldToFieldEndTimes);
            predictedVals=mod(predictedVals,360);
            behavTimeErrors=rad2ang(angdiff(ang2rad(actualVals),ang2rad(predictedVals)));



                    %meanErrorBehavTime=nanmean(abs(errors));
                    meanErrorBehavTime=circMeanDegNoMod360((behavTimeErrors));
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

                    title({sprintf('Theta phase vs Time in field, all fields (n=%d spikes)',length(allFieldSpikePhases)),sprintf('circR=%.2f, slope=%.2f deg/field, mean circ error = %.2f deg',RBehavTime,mBehavTime,meanErrorBehavTime)})
                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %error plotting
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                     subplot(2,2,3)

                    if(plotConditionalProb)
                        getJointDistrGivenX(allFieldToFieldEndTimes,behavTimeErrors,timeFieldEdges,phaseErrEdges,fTimeH);
                    else
                        getJointDistr(allFieldToFieldEndTimes,behavTimeErrors,timeFieldEdges,phaseErrEdges,fTimeH);
                    end

                     xlabel('time in field')
                    ylabel('Linear fit error (deg)')
                        hline(0,'k--',3)
                      box off
                     axis square
                      caxis(cProbLims)

                    title({sprintf('Linear fit error vs Time in field, all fields (n=%d spikes)',length(allFieldSpikePhases))})

                      subplot(2,2,4)

                    if(plotConditionalProb)
                        getJointDistrGivenX(allFieldToFieldEndLogTimes,logTimeErrors,timeFieldEdges,phaseErrEdges,fTimeH);
                    else
                        getJointDistr(allFieldToFieldEndLogTimes,logTimeErrors,timeFieldEdges,phaseErrEdges,fTimeH);
                    end


                     xlabel('log(time in field)')
                    ylabel('Linear fit error (deg)')
                     caxis(cProbLims)
                      box off
                     axis square

                     hline(0,'k--',3)

                    title({sprintf('Linear fit vs log(time in field), all fields (n=%d spikes)',length(allFieldSpikePhases))})
               end



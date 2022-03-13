function [slopeVsTime,offsetVsTime,RvsTime,pVsTime,newTimeAxis] = getLocalCorrelationStats(timeSeries,oldTimeAxis,movingWindSize)
        plotExamples=0;

	dt=oldTimeAxis(2)-oldTimeAxis(1);

	if(plotExamples)
		figure	
	end
        %plotExamples=1;

        %localWindSize=6;
        %localWindSize=12;
        localWindSize=round(movingWindSize/dt);
	%stepSize=round(localWindSize/3);
	stepSize=oldTimeAxis(2)-oldTimeAxis(1);
        %localWindSize=20;

	
        slopeVsTime=NaN(size(1:stepSize:length(timeSeries)));
        offsetVsTime=NaN(size(1:stepSize:length(timeSeries)));
        RvsTime=NaN(size(1:stepSize:length(timeSeries)));
        pVsTime=NaN(size(1:stepSize:length(timeSeries)));

	newTimeAxis=NaN(size(1:stepSize:length(timeSeries)));

	stepCount=0;
        for i=1:stepSize:length(timeSeries)
		stepCount=stepCount+1;
		newTimeAxis(stepCount)=oldTimeAxis(i);

		startIdx=max(1,i-round(localWindSize/2));
                endIdx=min(length(timeSeries),i+round(localWindSize/2));

                startTime=startIdx*dt;
                endTime=endIdx*dt;

                localValues=timeSeries(startIdx:endIdx);

                %if(length(localRate)>5)
                if(sum(~isnan(localValues))>5)
                        [m,b,R,p]=getLinearFitWithP(startIdx:endIdx,localValues);

                        slopeVsTime(stepCount)=m;
                        offsetVsTime(stepCount)=b;
                        RvsTime(stepCount)=R;
                        pVsTime(stepCount)=p;

                        if(plotExamples && sum(~isnan(localValues))==length(localValues))
                                plot(startIdx:endIdx,localValues,'ko','MarkerSize',10)

                                hold on
                                if(m<0)
                                        plot(xlim,m*xlim+b,'r--','LineWidth',5)
                                else
                                        plot(xlim,m*xlim+b,'g--','LineWidth',5)
                                end
                                xlabel('Time')
                                ylabel('Amp')

                                timeRangeStr=sprintf('%.2f_to_%.2f_sec,_R=%.2f,_b=%.2f',startIdx*stepSize,endIdx*stepSize,R,b);
				hold off
                                title(removeUnderscores(timeRangeStr))
				drawnow
                                %setFigFontTo(16)
                                %saveas(gcf,sprintf('%s_LocalLFPslope.tif',removeCommas(timeRangeStr)))
                        end
                end
        end

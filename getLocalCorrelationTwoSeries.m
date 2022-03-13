function [slopeVsTime,offsetVsTime,RvsTime,pVsTime,newTimeAxis] = getLocalCorrelationTwoSeries(timeSeries1,timeSeries2,oldTimeAxis,movingWindSize)
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
	stepSize=dt;
        %localWindSize=20;

	
        %slopeVsTime=NaN(size(1:stepSize:length(timeSeries1)));
        %offsetVsTime=NaN(size(1:stepSize:length(timeSeries1)));
        %RvsTime=NaN(size(1:stepSize:length(timeSeries1)));
        %pVsTime=NaN(size(1:stepSize:length(timeSeries1)));

	%newTimeAxis=NaN(size(1:stepSize:length(timeSeries1)));
    
    slopeVsTime=NaN(size(1:length(timeSeries1)));
        offsetVsTime=NaN(size(1:length(timeSeries1)));
        RvsTime=NaN(size(1:length(timeSeries1)));
        pVsTime=NaN(size(1:length(timeSeries1)));

	newTimeAxis=NaN(size(1:length(timeSeries1)));

	stepCount=0;
        %for i=1:stepSize:length(timeSeries1)
        for i=1:length(timeSeries1)
		stepCount=stepCount+1;
		newTimeAxis(stepCount)=oldTimeAxis(i);

		startIdx=max(1,i-round(localWindSize/2));
                endIdx=min(length(timeSeries1),i+round(localWindSize/2));

                startTime=startIdx*dt;
                endTime=endIdx*dt;

                localValues1=timeSeries1(startIdx:endIdx);
                localValues2=timeSeries2(startIdx:endIdx);

                %if(length(localRate)>5)
                if(sum(~isnan(localValues1))>0 && sum(~isnan(localValues2))>0)
                        [m,b,R,p]=getLinearFitWithP(localValues1,localValues2);

                        slopeVsTime(stepCount)=m;
                        offsetVsTime(stepCount)=b;
                        RvsTime(stepCount)=R;
                        pVsTime(stepCount)=p;

                        if(plotExamples && sum(~isnan(localValues1))==length(localValues1))
                                plot(localValues1,localValues2,'ko','MarkerSize',10)

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

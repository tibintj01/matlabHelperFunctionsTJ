function [estAccPerTimeStep,estJerkPerTimeStep] = getLocalAcc(timeAxis,cycleStartTimes,speedAlongTrackPerTime,gaussTimeAxis,localGaussSpikeRatePerTime)
    %use both spatial coordinates - care about motion other than along
    %track because it may modulate theta!!
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %load relevant variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %timeAxis=originalData.linear_track{1}.nonlinearFrames(:,1);
    %{
      timeAxis=originalData.frames(:,1);
    approxTimeStep=median(diff(timeAxis));
    taskStartTime=originalData.events(1,1);
    taskEndTime=originalData.events(2,1);
    withinTaskFrameIdxes=timeAxis<=taskEndTime & timeAxis>=taskStartTime;

    speedPerTimeStep=originalData.frames(withinTaskFrameIdxes,5);
     %speedPerTimeStep=originalData.linear_track{1}.nonlinearFrames(withinTaskFrameIdxes,5);
    speedPerTimeStep(speedPerTimeStep>200)=NaN; %spuriouss
    speedPerTimeStep=speedPerTimeStep(:);
    %}
    approxTimeStep=median(diff(timeAxis));
    speedPerTimeStep=speedAlongTrackPerTime;
    
    
    %create interpolated variables
    %speedInterpFactor=6;
    %speedInterpFactor=2;
   %interpolatedTimeAxis=interp1(1:length(timeAxis),timeAxis,linspace(1,length(timeAxis),speedInterpFactor*length(timeAxis)));
    %splineInterpSpeedPerTimeStep=interp1(timeAxis,speedPerTimeStep,interpolatedTimeAxis,'spline');
    
    %splineInterpSpeedPerTimeStep(splineInterpSpeedPerTimeStep<0)=NaN;
    
    
    %supportWindLength=5;
    polynomialFitOrder=2; %fit parabola-like curves
     supportWindLength=3;
      %polynomialFitOrder=4; %fit parabola-like curves
     %supportWindLength=5;
     
       %polynomialFitOrder=1; %long timescale acc and jerk
     %supportWindLength=0.1/approxTimeStep;
    %dtSpeed=[NaN; diff(timeAxis(:))];
        %smoothed acceleration ignorning steps
        
        polynomialFitOrder=1; %fit lines over acc time width
        accTimeWidth=0.5; %sec
        accTimeWidth=1; %sec
        accTimeWidth=0.2; %sec
        
        %accelerating approximated as a line or parabola over 0.5 seconds,
        %not interested in higher order fluctuations (seem to not be
        %predictive of phase precession slope)
        %polynomialFitOrder=2; 
        %accTimeWidth=0.5; %sec
        polynomialFitOrder=1; 
        accTimeWidth=0.3; %se
         %polynomialFitOrder=1; 
        %accTimeWidth=0.4; %sec
         %accTimeWidth=0.25; %sec
        %accTimeWidth=0.225; %sec
        
     supportWindLength=ceil(accTimeWidth/approxTimeStep);

     %assumes starts at time zero and has constant dt
    estAccPerTimeStep = movingslope(speedPerTimeStep,supportWindLength,polynomialFitOrder);
    estJerkPerTimeStep = movingslope(estAccPerTimeStep,supportWindLength,polynomialFitOrder);
    
    %accTimeAxis=timeAxis;
    %{
    figure;
    plot(timeAxis,speedPerTimeStep,'b-')
    hold on;
    plot(timeAxis,estAccPerTimeStep,'k-')
    disp('')
    %}
    
    %{
    accTimeAxis=timeAxis-min(timeAxis);%why the shift...?
    startIdx=1;
    while(accTimeAxis(startIdx)<min(timeAxis))
        if(accTimeAxis(startIdx)>min(timeAxis))
            break
        end
        startIdx=startIdx+1;
    end
    estAccPerTimeStep=estAccPerTimeStep(:);
    
    estAccPerTimeStep=[estAccPerTimeStep(startIdx:end);NaN(startIdx-1,1)];
    %}
    
    
    %{
    
    %timeWindowSmooth=0.1; %sec
    %halfWindIdx=round(timeWindowSmooth/approxTimeStep/2);
    
    %smoothedAccPerTimeStep=estAccPerTimeStep;
    
    estJerkPerTimeStep = movingslope(estAccPerTimeStep,supportWindLength,polynomialFitOrder);
    
   
    %smoothedJerkPerTimeStep=smoothdata(estJerkPerTimeStep,'movmean',halfWindIdx);
    smoothedJerkPerTimeStep=estJerkPerTimeStep;
    
    %splineInterpolatedJerk=interp1(timeAxis,estJerkPerTimeStep,interpolatedTimeAxis,'spline');
    %linearInterpolatedJerk=interp1(timeAxis,estJerkPerTimeStep,interpolatedTimeAxis,'linear');
    linearInterpolatedJerk=interp1(timeAxis,smoothedJerkPerTimeStep,interpolatedTimeAxis,'linear');

    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %compute moving slopes
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %plot relevant variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure
        yyaxis left
    plot(timeAxis,speedPerTimeStep,'b-','LineWidth',1); hold on; plot(timeAxis,speedPerTimeStep,'bo');
    
    %plot(interpolatedTimeAxis,splineInterpSpeedPerTimeStep,'k-'); hold on; plot(interpolatedTimeAxis,splineInterpSpeedPerTimeStep,'ko');
       ylim([0 70])
    
    yyaxis right
    %plot(timeAxis,estAccPerTimeStep,'k-','LineWidth',2); hold on; plot(timeAxis,estAccPerTimeStep,'ko');
        %plot(timeAxis,nancumsum(estAccPerTimeStep),'k-','LineWidth',2); hold on; plot(timeAxis,nancumsum(estAccPerTimeStep),'ko');
           %plot(timeAxis,nancumsum(nancumsum(estJerkPerTimeStep)),'k-','LineWidth',2); hold on; plot(timeAxis,nancumsum(nancumsum(estJerkPerTimeStep)),'ko');

        %ylim([0 50])
    %plot(timeAxis,estJerkPerTimeStep,'r-','LineWidth',3); hold on; plot(timeAxis,estJerkPerTimeStep,'ro');
        for cycleNum=1:length(cycleStartTimes)
              if(cycleStartTimes(cycleNum)>500)
                  %break
              end
              plot([cycleStartTimes(cycleNum) cycleStartTimes(cycleNum)],[min(estAccPerTimeStep) max(estAccPerTimeStep)],'k--','LineWidth',4)
            hold on
        end
        plot(gaussTimeAxis,localGaussSpikeRatePerTime,'r-','LineWidth',3)
          plot(gaussTimeAxis,localGaussSpikeRatePerTime,'ro')
         %{
            [pks,locs]=findpeaks(abs(estJerkPerTimeStep));
         
         for ji=1:length(locs)
             if(speedPerTimeStep(locs(ji))<10)
                 continue
             end
             peakTime=timeAxis(locs(ji));
             if(peakTime>500)
                 break
             end
             plot([peakTime peakTime],ylim,'m--','LineWidth',4)
         end
         
         %}
     %plot(interpolatedTimeAxis,abs(linearInterpolatedJerk),'m-','LineWidth',1); 
     plot(timeAxis,estAccPerTimeStep,'k-','LineWidth',2)
     %plot(interpolatedTimeAxis,(linearInterpolatedJerk)*10,'m-','LineWidth',2); 
     %dhold on; plot(interpolatedTimeAxis,abs(linearInterpolatedJerk),'mo');
    %plot(interpolatedTimeAxis,splineInterpolatedAcc,'r-'); hold on; plot(interpolatedTimeAxis,splineInterpolatedAcc,'ro');
    plot(xlim,[0 0],'k--')
    
    disp('')
    %}
    



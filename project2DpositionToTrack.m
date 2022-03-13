function [] = project2DpositionToTrack(timeAxis,posXYPerTimeStep)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %2D position over time!
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    approxTimeStep=median(diff(timeAxis));

    %posXPerTimeStep=scaledata(posXPerTimeStep,0,1)*trackLengthCm/(sqrt(2));
    %posYPerTimeStep=scaledata(posYPerTimeStep,0,1)*trackLengthCm/(sqrt(2));

    velocityXYPerTimeStep=NaN(2,length(posXPerTimeStep));


    %for ti=1:(length(timeAxis)-1)
    for ti=2:(length(timeAxis))
         %currTimeStep=timeAxis(ti+1)-timeAxis(ti);
         velocityXYPerTimeStep(1,ti)=(posXPerTimeStep(ti)-posXPerTimeStep(ti-1))/currTimeStep;
        velocityXYPerTimeStep(2,ti)=(posYPerTimeStep(ti)-posYPerTimeStep(ti-1))/currTimeStep;

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %project position and velocity onto vector along track [1 ;1]/sqrt(2)
    %since track seems oriented at 45 degree angle given x and y
    %coordinates
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [speedAlongTrackPerTime] = projectDataOntoVectorSpan(velocityXYPerTimeStep,[1 1]/sqrt(2));
    [posAlongTrackPerTime] = projectDataOntoVectorSpan(posXYPerTimeStep,[1 1]/sqrt(2));

    posAlongTrackPerTime=scaledata(posAlongTrackPerTime,0,1)*trackLengthCm;

    %speedAlongTrackPerTime=scaledata(abs(speedAlongTrackPerTime),0,max(speedPerTimeStep));
    speedAlongTrackPerTime=abs(speedAlongTrackPerTime)/max(abs(speedAlongTrackPerTime))*max(speedPerTimeStep);%convert to same max

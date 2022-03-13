function [placeCellDistPerTime,placeCellTimeInLapPerTime,placeCellRatePerTime,comDist,comTime] = addSpeedChangeToDataStruct(filePath)


   
    
    %zeroFieldToCOM=1;
    zeroFieldToCOM=0;
  
    data=load(filePath);
    
    data
    allSpikeTimes=data.spikePerCycleInfo.allSpikeTimes;
    manualPrecessionStartPhaseDeg=data.manualPrecessionStartPhaseDeg;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %choose direction
    %From code generating file name with "dir"
    %        if(di==1)
    %            dirFlag=-1;
    %        else
    %            dirFlag=1;
    %        end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dirFlag=data.spaceTimePhaseInfo.dirFlag;
    if(dirFlag==-1)
        manualSelectedDI=1;
    elseif(dirFlag==1)
        manualSelectedDI=2;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %SET PLACE FIELD BOUNDS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    manualFieldStartCm=data.manualFieldStartCm;
    manualFieldEndCm=data.manualFieldEndCm;
    
    %manualFieldStartCm=data.manualFieldStartCmAdjusted;
    %manualFieldEndCm=data.manualFieldEndCmAdjusted;
    
    minSpeedForChange=3; %cm/s
    fileRootName=getFileNameFromPath(data.spikePerCycleInfo.filePath);
    fileRootName=fileRootName(1:(end-4));
    originalDataDir='/Users/tibinjohn/Downloads/processed_data/';
    %{
    if(~contains(data.spikePerCycleInfo.filePath,'.mat'))
        originalDataPath=[data.spikePerCycleInfo.filePath '.mat'];
    else
        originalDataPath=data.spikePerCycleInfo.filePath;
    end
    %}
    
    matchingFilePaths=getFilePathsRegex(originalDataDir,[fileRootName '*.mat']);
    
       originalData=load(matchingFilePaths{1});
    %originalData=load(originalDataPath);
    
 
     
    cellNum=data.spikePerCycleInfo.cellNum;
    
    avgWaveform=originalData.avgwave{cellNum};
    %figure; plot(avgWaveform')
    
    trackLengthCm=originalData.maze_size_cm(1); %linear track length
    timeAxis=originalData.frames(:,1);
    taskStartTime=originalData.events(1,1);
    taskEndTime=originalData.events(2,1);
    withinTaskFrameIdxes=timeAxis<=taskEndTime & timeAxis>=taskStartTime;
    
    posPerTimeStep=originalData.frames(withinTaskFrameIdxes,2);
    posPerTimeStep=scaledata(posPerTimeStep,0,1);
    %speedPerTimeStep=originalData.frames(withinTaskFrameIdxes,5);
    
     speedPerTimeStep=originalData.frames(withinTaskFrameIdxes,5);
     noiseSpeedThresh=200;
    speedPerTimeStep(speedPerTimeStep>noiseSpeedThresh)=NaN;
   %{
    figure
    plot(posPerTimeStep,speedPerTimeStep,'k')
    hold on
    plot(posPerTimeStep,speedPerTimeStep,'ko')
    
    %}
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %2D position over time!
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    timeAxis=originalData.linear_track{1}.nonlinearFrames(:,1);
    withinTaskFrameIdxes=timeAxis<=taskEndTime & timeAxis>=taskStartTime;
    timeAxis=timeAxis(withinTaskFrameIdxes);
    approxTimeStep=median(diff(timeAxis));
    posXPerTimeStep=originalData.linear_track{1}.nonlinearFrames(withinTaskFrameIdxes,2);
    posYPerTimeStep=originalData.linear_track{1}.nonlinearFrames(withinTaskFrameIdxes,3);
    
    posXPerTimeStep=scaledata(posXPerTimeStep,0,1)*trackLengthCm/(sqrt(2));
    posYPerTimeStep=scaledata(posYPerTimeStep,0,1)*trackLengthCm/(sqrt(2));
    
    velocityXYPerTimeStep=NaN(2,length(posXPerTimeStep));
    posXYPerTimeStep=[posXPerTimeStep(:)';posYPerTimeStep(:)'];
    
   
    %for ti=1:(length(timeAxis)-1)
    for ti=2:(length(timeAxis))
         %currTimeStep=timeAxis(ti+1)-timeAxis(ti);
         currTimeStep=1/30; %sec
        %velocityXYPerTimeStep(1,ti)=(posXPerTimeStep(ti+1)-posXPerTimeStep(ti))/currTimeStep;
        %velocityXYPerTimeStep(2,ti)=(posYPerTimeStep(ti+1)-posYPerTimeStep(ti))/currTimeStep;
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %local spike rate (~100msec time window) per time - must be small to
    %not interfere with time dimension of spacetime
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    localTimeWindow=0.15; %sec
    
    localSpikeRatePerTime=getLocalWindowEventRate(allSpikeTimes,timeAxis,localTimeWindow);
    
    %spikeBoolean=rasterize(allSpikeTimes,1/approxTimeStep);
    
    %[gaussTimeAxis,localGaussSpikeRatePerTime]=localGaussRateSpikeTrain(allSpikeTimes, 0.05,min(timeAxis),max(timeAxis));
        
    [gaussTimeAxis,localGaussSpikeRatePerTime]=localGaussRateSpikeTrain(allSpikeTimes, 0.1,0,max(timeAxis));
    %[gaussTimeAxis,localGaussSpikeRatePerTime]=localGaussRateSpikeTrain(allSpikeTimes, 0.05,0,max(timeAxis));
        %[gaussTimeAxis,localGaussSpikeRatePerTime]=localGaussRateSpikeTrain(allSpikeTimes, 0.1,min(timeAxis),max(timeAxis));

    
    
    %localGaussSpikeRatePerTime=interp1(gaussTimeAxis,localGaussSpikeRatePerTime,timeAxis);
    
    %{
    figure;
      yyaxis left
    plot(timeAxis,localSpikeRatePerTime,'b')
  
    yyaxis right
    plot(timeAxis,localGaussSpikeRatePerTime,'k')
   %}
    allSpikeTimesAndPhases=getAllSpikePhasesAndTimes(originalData,data);
    
    localTimeWindowPhase=0.1;%like theta cycle, 100msec window to detect phase
    %localTimeWindowPhase=0.12;%like theta cycle, 120msec window to detect phase
     %localTimeWindowPhase=0.05; %like gamma cycle, 30msec window to detect phase
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %get local average phase based on average phase in local fixed window
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %localSpikePhasePerTime=getLocalPhasePerTime(allSpikeTimesAndPhases,timeAxis,localTimeWindowPhase);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %get local average phase based on average phase in local theta cycle
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cycleStartTimes=data.spikePerCycleInfo.cycleStartTimes;
    cycleAmpsLFPZ=data.spikePerCycleInfo.thetaAmpZPerCycle;
    cycleDurationsSec=data.spikePerCycleInfo.durationPerCycleSec;
    
    %calculateJerkAccurately(originalData,cycleStartTimes,gaussTimeAxis,localGaussSpikeRatePerTime);
     [localAccPerTime,localJerkPerTime]=getLocalAcc(timeAxis,cycleStartTimes,speedAlongTrackPerTime,gaussTimeAxis,localGaussSpikeRatePerTime);
     
    

    [localCycleMeanSpikePhasePerTime, localCycleMeanSpikeTimePerTime]=getLocalCycleMeanPhasePerTime(allSpikeTimesAndPhases,timeAxis,cycleStartTimes);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %set local average phase based on average phase or peak of gauss rate!
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gaussBasedPhasePerCycle=1;
    if(gaussBasedPhasePerCycle)
         localSpikePhasePerTime=getPeakPhasePerThetaCyclePerTime(gaussTimeAxis,localGaussSpikeRatePerTime,cycleStartTimes,timeAxis);
    else
        localSpikePhasePerTime=localCycleMeanSpikePhasePerTime;
    end
    
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %MANUAL OFFSET STANDARDIZATION
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    localSpikePhasePerTime= mod(localSpikePhasePerTime-(manualPrecessionStartPhaseDeg-360),360);
    
    %{
     figure;
      yyaxis left
    plot(timeAxis,localAccPerTime,'k')
    hold on
     plot(gaussTimeAxis,localGaussSpikeRatePerTime,'r-')
  
    yyaxis right
    %plot(timeAxis,localSpikePhasePerTime/360,'m')
    %hold on
   
    hold on
    plot(timeAxis,abs(speedAlongTrackPerTime),'b-')
    %plot(timeAxis,speedPerTimeStep,'m-')
    %}
    
  
    
    localSpikeTimePerTime=localCycleMeanSpikeTimePerTime;
    
    localSpikePhasePerTime=localSpikePhasePerTime(:);
    localSpikeTimePerTime=localSpikeTimePerTime(:);
    
    %localSpikePhaseDiffPerTime=[diff(localSpikePhasePerTime); NaN];
    localSpikePhaseDiffPerTime=[angdiffDeg(localSpikePhasePerTime); NaN];
    %stepHalfWind=3;
    %smoothedSpikePhasePerTime=circSmoothDeg(localSpikePhasePerTime,stepHalfWind);
    
    %localSpikePhaseDiffPerTime
    
    [localCycleAmpPerTime,localCycleDurPerTime,localCycleNumPerTime,localCycleStartTimePerTime]=getLocalCycleMeanAmpAndDurPerTime(timeAxis,cycleStartTimes,cycleAmpsLFPZ,cycleDurationsSec);
    
    goodThetaIdxes=localCycleAmpPerTime>2;
     %goodThetaIdxes=localCycleAmpPerTime>2.5;
    
    localCycleMeanSpikePhasePerTime(~goodThetaIdxes)=NaN;
    
    %figure; plot(timeAxis,localSpikePhasePerTime)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %use scale factor that gives best match with provided speed
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %maxVFwdDir=
    %maxVBkwdDir=
    
    providedLapStartTimes=[originalData.linear_track{1}.lapinfo.laps.start_ts];
    providedLapStartPosition=[originalData.linear_track{1}.lapinfo.laps.pos];
    
    absSpeedAlongTrackPerTime=abs(speedAlongTrackPerTime);
    absSpeedAlongTrackPerTime=scaledata(absSpeedAlongTrackPerTime,0,max(speedPerTimeStep));

    
    diffAbsSpeedPerTime=[diff(absSpeedAlongTrackPerTime(:)); NaN];
    %0 to 120 increasing one way, 120 to 0 decreasing otherway (don't
    %distinguish--> absolute value)
    diffTrackPosPerTime=abs([diff(posAlongTrackPerTime(:)); NaN]);
    
    speedChangePerPosChangePerTime=diffAbsSpeedPerTime./diffTrackPosPerTime;
    speedChangePerTimePerTime=diffAbsSpeedPerTime(:)./[diff(timeAxis(:));NaN];
    
    speedChangePerPosChangePerTime=localAccPerTime(:);
   
    %if(manualSelectedDI==1)
        %speedChangePerPosChangePerTime=-speedChangePerPosChangePerTime;
        %speedChangePerTimePerTime=-speedChangePerTimePerTime;
    %end
    
    maxSpeedChangeMagnitude=5;
    maxSpeedChangeMagnitude=10;
    speedChangePerPosChangePerTime(abs(speedChangePerPosChangePerTime)>maxSpeedChangeMagnitude)=NaN;
    %{
    interpolatedTimeAxis=interp1(1:length(timeAxis),timeAxis,linspace(1,length(timeAxis),3*length(timeAxis)));
      linearInterpSpeedChangePerPosPerTime=interp1(timeAxis,speedChangePerPosChangePerTime,interpolatedTimeAxis,'linear');
    splineInterpSpeedChangePerPosPerTime=interp1(timeAxis,speedChangePerPosChangePerTime,interpolatedTimeAxis,'spline');
    linearInterpDistChangePerTime=interp1(timeAxis,diffTrackPosPerTime,interpolatedTimeAxis,'linear');
     
    splineInterpSpeedChangePerPosPerTime(isnan(linearInterpSpeedChangePerPosPerTime))=NaN;
    splineInterpSpeedChangePerPosPerTime(abs(splineInterpSpeedChangePerPosPerTime)>maxSpeedChangeMagnitude)=NaN;
    %speedChangePerTimePerTime=scaledata(-5,5,speedChangePerTimePerTime);
    
    %diffSpeedChangePerPosPerTime=[diff(splineInterpSpeedChangePerPosPerTime(:)); NaN];
    %diffDistChangePerTime=[diff(linearInterpDistChangePerTime(:));NaN];
    %diffInterpolatedTime=[diff(interpolatedTimeAxis(:));NaN];
    %estJerkPerPosPerTime=diffSpeedChangePerPosPerTime./diffDistChangePerTime; %percent acceleration change???
     %estJerkPerPosPerTime=diffSpeedChangePerPosPerTime./diffInterpolatedTime; %percent acceleration change???
  
    %estJerkPerPosPerTime(abs(estJerkPerPosPerTime)>20)=NaN;
    
    
    dispIdxesSpeedChange=1:length(timeAxis);
        figure;
        yyaxis left; plot(timeAxis(dispIdxesSpeedChange),absSpeedAlongTrackPerTime(dispIdxesSpeedChange),'b-'); hold on; plot(timeAxis(dispIdxesSpeedChange),absSpeedAlongTrackPerTime(dispIdxesSpeedChange),'bo');
       
                hold on
        ylim([0 max(absSpeedAlongTrackPerTime)])
          yyaxis right;    
                  %plot(timeAxis(dispIdxesSpeedChange),speedChangePerPosChangePerTime(dispIdxesSpeedChange),'r-'); hold on; plot(timeAxis(dispIdxesSpeedChange),speedChangePerPosChangePerTime(dispIdxesSpeedChange),'ro');
                  plot(interpolatedTimeAxis,splineInterpSpeedChangePerPosPerTime,'m-'); hold on; plot(interpolatedTimeAxis,splineInterpSpeedChangePerPosPerTime,'mo');
                    %plot(interpolatedTimeAxis,estJerkPerPosPerTime,'r-'); hold on;plot(interpolatedTimeAxis,estJerkPerPosPerTime,'ro'); 
                    
          for cycleNum=1:length(cycleStartTimes)
              if(cycleStartTimes(cycleNum)>timeAxis(dispIdxesSpeedChange(end)))
                  continue
              end
              plot([cycleStartTimes(cycleNum) cycleStartTimes(cycleNum)],[min(speedChangePerPosChangePerTime) max(speedChangePerPosChangePerTime)],'k--','LineWidth',5)
          end
          plot(xlim,[0 0],'k--')
          
         plot( gaussTimeAxis,localGaussSpikeRatePerTime,'r-')
          
          %{
          for s=1:length(allSpikeTimes)
                 if(allSpikeTimes(s)>timeAxis(dispIdxesSpeedChange(end)))
                  continue
                end
                     plot([allSpikeTimes(s) allSpikeTimes(s)],[min(speedChangePerPosChangePerTime) max(speedChangePerPosChangePerTime)],'r-','LineWidth',2)
          end
         %}
          
          %plot moving window spike rate
          
          
        %plot(timeAxis(dispIdxesSpeedChange),speedChangePerTimePerTime(dispIdxesSpeedChange),'m-');hold on; plot(timeAxis(dispIdxesSpeedChange),speedChangePerTimePerTime(dispIdxesSpeedChange),'mo')
        %plot(timeAxis(dispIdxesSpeedChange),20*diffTrackPosPerTime(dispIdxesSpeedChange),'k-'); hold on; plot(timeAxis(dispIdxesSpeedChange),20*diffTrackPosPerTime(dispIdxesSpeedChange),'ko'); 
        
        for li=1:length(providedLapStartTimes)
            xlim([providedLapStartTimes(li) providedLapStartTimes(li)+3])
            figure(gcf)
        end
        
        %cell LEM3206_S2019072221142_cellNum10_dir1
        xlim([42.3 43.3])
        xlim([46.75 48.25])
        xlim([50.8 52])
        xlim([54.5 56])
        xlim([60.75 62.1])
        xlim([67.5 68.5])
        %first cell (LEM3116 cell 20 dir 1)
        xlim([240 241.25])
        xlim([258.5 259.5])
        xlim([280.75 281.75])
             xlim([298.5 299])
        xlim([317.25 318])
        xlim([333.1 334])
        xlim([346 346.8])
        xlim([366.2 366.85])
        xlim([382.7 383.4])
        xlim([457.75 458.75 ])
        xlim([515 515.5])
         xlim([524 524.7]) 
         xlim([578 578.75])%negativity of deltaSpeed determining phase jump (frequency increase)
         xlim([639 640])
         xlim([692.25 692.9])
         xlim([732.75 734.2])
         xlim([769.8 770.45])
         xlim([839.4 840.5])
         xlim([948.8 949.3])
        
   
        
        
        xlim([169 171])
        xlim([147.5 150])
        
        
      %ylim([-5 5])
    
    %}
    speedChangePerPosChangePerTime(absSpeedAlongTrackPerTime<10)=NaN;
  
    

    %maxSpeedChangePerCm=10;
    %speedChangePerPosChangePerTime(abs(speedChangePerPosChangePerTime)>maxSpeedChangePerCm)=NaN;
    %speedChangePerPosChangePerTime(absSpeedAlongTrackPerTime<minSpeedForChange)=NaN;
    
    speedAlongTrackPerTime=scaledata((speedAlongTrackPerTime),-max(speedPerTimeStep),max(speedPerTimeStep));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %find true lap start times - start at 60cm of existing lap and go
    %true start: backwards in time until hits less than 5cm/s first time before 20cm to end
    %true end: forwards in time until hits less than 5cm/s first time after 20cm to end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    trueLapStartTimes=NaN(size(providedLapStartTimes));
    trueLapEndTimes=NaN(size(providedLapStartTimes));
    
    %trueLapStartPositions=NaN(size(providedLapStartTimes));
    %maxLapStartPosLow=20; %cm
    maxLapStartPosLow=15; %cm
    maxLapStartPosLow=30; %cm - eliminate more pauses to get truer lap start
    maxLapStartPosLow=5; %cm - to get early fields?
    maxLapStartPosHigh=trackLengthCm-maxLapStartPosLow;
    
    minSpeedLapStart=5;   %cm/s
    %minSpeedLapEnd=1;   %cm/s - must get very close to zero to turn around
    minSpeedLapEnd=5;   %cm/s - must get very close to zero to turn around
    
    maxLapDuration=10; %seconds (assuming such long pause during lap is rare/not strongly coded)
    maxLapDuration=5;
    maxTimeIdxDuringLap=ceil(maxLapDuration/approxTimeStep);
    
    isDuringTrueLap=false(size(timeAxis));
    trueLapNumPerTime=NaN(size(timeAxis));
    timeSinceLastTrueLapStartPerTime=NaN(size(timeAxis));
    lapDirPerTime=NaN(size(timeAxis));

    numLaps=length(providedLapStartTimes);
    trueLapDurations=NaN(numLaps,1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %find true lap start time for each lap
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for li=1:numLaps
        currLapStartPos=providedLapStartPosition(li);
        currLapStartTime=providedLapStartTimes(li);
        
        if(li<length(providedLapStartTimes))
            currLapEndTime=providedLapStartTimes(li+1);
        else
            currLapEndTime=max(timeAxis);
        end
        
        currLapIdxes=timeAxis>=currLapStartTime & timeAxis<currLapEndTime;
        startingTimeIdxOffset=length(find(timeAxis<=currLapStartTime))-1;
        
        currLapPositions=posAlongTrackPerTime(currLapIdxes);
        currLapSpeeds=speedAlongTrackPerTime(currLapIdxes);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %find true lap start time by starting at the middle and going back  
        %in time until speed drops below min speed for the first time
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if(currLapStartPos<trackLengthCm/2)
            [~,currSearchID]=min(abs(currLapPositions-maxLapStartPosLow));
        else
            [~,currSearchID]=min(abs(currLapPositions-maxLapStartPosHigh));
        end
        
        while(abs(currLapSpeeds(currSearchID))>minSpeedLapStart)
           if(currSearchID==1)
                break;
            end
            currSearchID=currSearchID-1;
        end
        correspondingStartIdx=currSearchID+startingTimeIdxOffset;
        trueLapStartTimes(li)=timeAxis(correspondingStartIdx); %accurate to within closest 30 msec interval
        
        trueLapStartPositions(li)=posAlongTrackPerTime(correspondingStartIdx);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %find true lap end time by starting at the middle and going forward  
        %in time until speed drops below min speed for the first time
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if(currLapStartPos<trackLengthCm/2)
            [~,currSearchID]=min(abs(currLapPositions-(maxLapStartPosHigh)));
        else
            [~,currSearchID]=min(abs(currLapPositions-(maxLapStartPosLow)));
        end
        
        while(abs(currLapSpeeds(currSearchID))>minSpeedLapEnd)
            if(currSearchID==length(currLapSpeeds))
                break;
            end
            currSearchID=currSearchID+1;

          
        end
        correspondingEndIdx=currSearchID+startingTimeIdxOffset;
        trueLapEndTimes(li)=timeAxis(correspondingEndIdx);

        isDuringTrueLap(correspondingStartIdx:correspondingEndIdx)=true;
        trueLapNumPerTime(correspondingStartIdx:correspondingEndIdx)=li;
      
        
        trueLapEndPositions(li)=posAlongTrackPerTime(correspondingEndIdx);
        
        trueLapDurations(li)=trueLapEndTimes(li)-trueLapStartTimes(li);
        numTimeStepsInLap=length(correspondingStartIdx:correspondingEndIdx);
        
        %if(providedLapStartPosition(li)<trackLengthCm/2)
        if(providedLapStartPosition(li)>trackLengthCm/2)
            currDirID=1;
        else
            currDirID=2;
        end
        
        lapDirPerTime(correspondingStartIdx:correspondingEndIdx)=repelem(currDirID,numTimeStepsInLap);
        %timeSinceLastTrueLapStartPerTime(correspondingStartIdx:correspondingEndIdx)=(approxTimeStep/2):approxTimeStep:(trueLapDurations(li)+approxTimeStep/2);
         timeSinceLastTrueLapStartPerTime(correspondingStartIdx:correspondingEndIdx)=linspace(0,trueLapDurations(li),numTimeStepsInLap);
    end
    
    dirIdxes=lapDirPerTime==manualSelectedDI;
    
    if(manualSelectedDI==1)
        manualFieldStartCm=trackLengthCm-manualFieldStartCm;
        manualFieldEndCm=trackLengthCm-manualFieldEndCm;
    end
    speedChangePerPosChangePerTime(~isDuringTrueLap)=NaN;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %find field entry/exit time for each lap
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fieldEntryTime=NaN(numLaps,1);
    fieldEntryPos=NaN(numLaps,1);
    fieldExitTime=NaN(numLaps,1);
    fieldExitPos=NaN(numLaps,1);
    timeSinceLastFieldEntry=NaN(size(timeAxis));
    distanceFieldFracPerTime=NaN(size(timeAxis));
    
     timeSinceLastFieldEntryNormLap=NaN(size(timeAxis));
    distanceFieldFracPerTimeNormLap=NaN(size(timeAxis));
    thisLapHasBackwardsMovementPerTime=false(size(timeAxis));
    %maxSpeedField=100; %cm/sec
    %maxSpeedField=150;
    %minSpeedField=10; %cm/sec
    maxSpeedField=200;
    %minSpeedField=5; %cm/sec
    minSpeedField=0; %cm/sec
    
    %maxSpeedField=200; %cm/sec
    %minSpeedField=5; %cm/sec
     %minSpeedField=20;%cm/sec
    
    minTimeInField=abs(manualFieldEndCm- manualFieldStartCm)/maxSpeedField;
    %maxFieldTime=2.5;
    maxFieldTime=abs(manualFieldEndCm- manualFieldStartCm)/minSpeedField;
    
    closenessTolerance=1; %cm
   
    peakSpikeRatePerLap=NaN(numLaps,1);
    fieldWidthPerLap=NaN(numLaps,1);
    timeInFieldPerLap=NaN(numLaps,1);
    averageStartPhasePerLap=NaN(numLaps,1);
    averageEndPhasePerLap=NaN(numLaps,1);
    speedAtStartPerLap=NaN(numLaps,1);
    speedAtEndPerLap=NaN(numLaps,1);
    phasesPerLap=[];
 
     absSpeedsDuringField=[];
     phaseDiffPerTimeDuringField=[];
     allElapsedTimesInFieldSec=[];
     allElapsedDistsInFieldCm=[];
     allTimesInField=[];
     allCycleDursInField=[];
     allPhasesInField=[];
     allLapNumsInField=[];
     allLocalSpikeTimesFromFieldEntry=[];
     allCycleNumsInField=[];
     allCycleStartTimesInField=[];
     speedChangesPerPosDuringField=[];
     jerkDuringField=[];
 
    throwoutDistThresh=abs(manualFieldEndCm- manualFieldStartCm)/5;
     %throwoutDistThresh=abs(manualFieldEndCm- manualFieldStartCm)/4;
      %throwoutDistThresh=Inf; %throw out none
      
      NaNSeparator=NaN(20,1);
    for li=1:numLaps
        currLapEndTime=trueLapEndTimes(li);
        currLapStartTime=trueLapStartTimes(li);
        
         currLapEndPos=trueLapEndPositions(li);
         currLapStartPos=trueLapStartPositions(li);
        
        currLapIdxes=timeAxis>=currLapStartTime & timeAxis<currLapEndTime;
        startingTimeIdxOffset=length(find(timeAxis<=currLapStartTime))-1;
        
        currLapPositions=posAlongTrackPerTime(currLapIdxes);
        
       
        currLapSpeeds=speedAlongTrackPerTime(currLapIdxes);
        
        currLapPositions(abs(currLapSpeeds)<2)=NaN;
        
        %[~,currSearchID]=min(abs(currLapPositions-(manualFieldStartCm)));
          closeIdxes=find(abs(currLapPositions-(manualFieldStartCm))<closenessTolerance);
          %if(currLapStartPos<trackLengthCm/2) %start is lower position than end
           %         currSearchID=min(closeIdxes); %latest index that is close enough
            % else
             %    currSearchID=max(closeIdxes); %earliest index that is close enough
          %end
             
                if(manualSelectedDI==2)
                       currSearchID=max(closeIdxes); %earliest index that is close enough
                   else
                    currSearchID=max(closeIdxes); %earliest index that is close enough
                   end
          
          if(isempty(currSearchID))
              [~,currSearchID]=min(abs(currLapPositions-(manualFieldStartCm)));
          end
          
   
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %find field entry time by starting at the end and going backward  
        %in time until position drops below start for the first time
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %{
        [~,currSearchID]=min(abs(currLapPositions-(currLapEndPos)));
            
        if(currLapStartPos<trackLengthCm/2) %start is lower position than end
            while(abs(currLapPositions(currSearchID))>manualFieldStartCm)
                if(currSearchID==1)
                    break;
                end
                currSearchID=currSearchID-1;       
             end
        else                                %start is higher position than end
             while(abs(currLapPositions(currSearchID))<manualFieldStartCm)
                if(currSearchID==1)
                    break;
                end
                currSearchID=currSearchID-1;       
             end
        end
        %}
        
        correspondingEntryIdx=currSearchID+startingTimeIdxOffset;
        
       
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %find field exit time by starting at the start and going forward  
        %in time until position rises above field end for the first time
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
          %[~,currSearchID]=min(abs(currLapPositions-(manualFieldEndCm)));
        closeIdxes=find(abs(currLapPositions-(manualFieldEndCm))<closenessTolerance);
             %if(currLapStartPos<trackLengthCm/2) %start is lower position than end
                   if(manualSelectedDI==2)
                       currSearchID=min(closeIdxes); %earliest index that is close enough
                   else
                    currSearchID=min(closeIdxes); %earliest index that is close enough
                   end
            % else
             %    currSearchID=min(closeIdxes); %earliest index that is close enough
             %end
        if(isempty(currSearchID))
            [~,currSearchID]=min(abs(currLapPositions-(manualFieldEndCm)));
        end
            
        %{
          [~,currSearchID]=min(abs(currLapPositions-(currLapStartPos)));
        if(currLapStartPos<trackLengthCm/2) %start is lower position than end
            while(abs(currLapPositions(currSearchID))<manualFieldEndCm)
                if(currSearchID==length(currLapPositions))
                    break;
                end
                currSearchID=currSearchID+1;       
             end
        else                                %start is higher position than end
             while(abs(currLapPositions(currSearchID))>manualFieldEndCm)
                if(currSearchID==length(currLapPositions))
                    break;
                end
                currSearchID=currSearchID+1;       
             end
        end
        %}
        
        correspondingExitIdx=currSearchID+startingTimeIdxOffset;
        
        if(correspondingExitIdx<correspondingEntryIdx)
            tempVar=correspondingEntryIdx;
            correspondingEntryIdx=correspondingExitIdx;
            correspondingExitIdx=tempVar;
            
        end
        
        timeInField=(correspondingExitIdx-correspondingEntryIdx+1)*approxTimeStep; %timestep is every 30msec
       
        
        numTimeIdxesInField=(correspondingExitIdx-correspondingEntryIdx+1);
        
        if(isempty(timeInField) || timeInField<minTimeInField || timeInField>maxFieldTime)
            continue
        end
        
        if(abs(posAlongTrackPerTime(correspondingEntryIdx) - manualFieldStartCm) >throwoutDistThresh)
            continue
        end
        
        if(abs(posAlongTrackPerTime(correspondingExitIdx) - manualFieldEndCm) >throwoutDistThresh)
            continue
        end
        
       
        
        if(~isempty(correspondingEntryIdx) && ~isempty(correspondingExitIdx))
             posSeqInCurrLap=posAlongTrackPerTime(correspondingEntryIdx:correspondingExitIdx);
             if((manualSelectedDI==2 && min(diff(posSeqInCurrLap)) < 0) || (manualSelectedDI==1 && max(diff(posSeqInCurrLap)) > 0))
                 thisLapHasBackwardsMovementPerTime(correspondingEntryIdx:correspondingExitIdx)=true;
                 continue %skip laps where the rodent goes the wrong way (different place map)
             end
            
            %timeSinceLastFieldEntry(correspondingEntryIdx:correspondingExitIdx)=linspace(0,timeInField,numTimeIdxesInField);
            timeSinceLastFieldEntryNormLap(correspondingEntryIdx:correspondingExitIdx)=linspace(0,timeInField,numTimeIdxesInField)/timeInField; %NORMALIZE TIME IN FIELD FOR CELL TO COMPARE FIELDS
            timeSinceLastFieldEntry(correspondingEntryIdx:correspondingExitIdx)=linspace(0,timeInField,numTimeIdxesInField); %NORMALIZE TIME IN FIELD FOR CELL TO COMPARE FIELDS
            
            fieldEntryTime(li)=timeAxis(correspondingEntryIdx);
            fieldEntryPos(li)=posAlongTrackPerTime(correspondingEntryIdx);

            fieldExitTime(li)=timeAxis(correspondingExitIdx);
            fieldExitPos(li)=posAlongTrackPerTime(correspondingExitIdx);  
            
            fieldWidth=abs(fieldExitPos(li)-fieldEntryPos(li));
            
            speedAtStartPerLap(li)=absSpeedAlongTrackPerTime(correspondingEntryIdx);
            speedAtEndPerLap(li)=absSpeedAlongTrackPerTime(correspondingExitIdx);
           
	     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             %NORMALIZE DIST, TIME IN CELL AVG FIELD TO COMPARE FIELDS
	     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            distanceFieldFracPerTimeNormLap(correspondingEntryIdx:correspondingExitIdx)=abs(posAlongTrackPerTime(correspondingEntryIdx:correspondingExitIdx)-manualFieldStartCm)/fieldWidth;
            distanceFieldFracPerTime(correspondingEntryIdx:correspondingExitIdx)=abs(posAlongTrackPerTime(correspondingEntryIdx:correspondingExitIdx)-manualFieldStartCm);
            %distanceFieldFracPerTime(correspondingEntryIdx:correspondingExitIdx)=abs(posAlongTrackPerTime(correspondingEntryIdx:correspondingExitIdx)-manualFieldStartCm);
           

                
         currLapDir=median(lapDirPerTime(correspondingEntryIdx:correspondingExitIdx));
         
            if(currLapDir==manualSelectedDI)
                ratesInLap=localSpikeRatePerTime(correspondingEntryIdx:correspondingExitIdx);
                peakSpikeRatePerLap(li)=nanmax(ratesInLap);
                
                 fieldWidthPerLap(li)=fieldWidth;
                  timeInFieldPerLap(li)=timeInField;
                  earlyPhasesInLap=localSpikePhasePerTime(correspondingEntryIdx:(correspondingEntryIdx+10));
                   latePhasesInLap=localSpikePhasePerTime((correspondingExitIdx-10):(correspondingExitIdx));
                   
                   
                   currLapPhases=localSpikePhasePerTime(correspondingEntryIdx:correspondingExitIdx);
                  
                   %currLapPhases= mod(currLapPhases-(manualPrecessionStartPhaseDeg-360),360);
                   phasesPerLap{li}=currLapPhases;
                     
        %earlyPhasesInLap= mod(earlyPhasesInLap-(manualPrecessionStartPhaseDeg-360),360);
        %latePhasesInLap= mod(latePhasesInLap-(manualPrecessionStartPhaseDeg-360),360);
        
                  averageStartPhasePerLap(li)=nanmean(earlyPhasesInLap);
                  averageEndPhasePerLap(li)=nanmean(latePhasesInLap);
                  
                  currFieldPhaseDiffs=localSpikePhaseDiffPerTime(correspondingEntryIdx:correspondingExitIdx);
                  currLapNums=trueLapNumPerTime(correspondingEntryIdx:correspondingExitIdx);
                 
                  
                  currFieldSpeeds=absSpeedAlongTrackPerTime(correspondingEntryIdx:correspondingExitIdx);
                  currFieldElapsedTimes=timeSinceLastFieldEntry(correspondingEntryIdx:correspondingExitIdx);
                  currFieldElapsedDists=distanceFieldFracPerTime(correspondingEntryIdx:correspondingExitIdx);
                  currFieldTimes=timeAxis(correspondingEntryIdx:correspondingExitIdx);
                  currCycleDurs=localCycleDurPerTime(correspondingEntryIdx:correspondingExitIdx);
                  currFieldCycleNums=localCycleNumPerTime(correspondingEntryIdx:correspondingExitIdx);
                  currFieldCycleStartTimes=localCycleStartTimePerTime(correspondingEntryIdx:correspondingExitIdx);

                  currSpikeTimesFromFieldEntry=localSpikeTimePerTime(correspondingEntryIdx:correspondingExitIdx)-timeAxis(correspondingEntryIdx);
                  currSpeedChanges=speedChangePerPosChangePerTime(correspondingEntryIdx:correspondingExitIdx);
                  currJerks=localJerkPerTime(correspondingEntryIdx:correspondingExitIdx);
                  
                  absSpeedsDuringField=[absSpeedsDuringField(:);currFieldSpeeds(:);NaNSeparator];
                  
                  speedChangesPerPosDuringField=[speedChangesPerPosDuringField(:); currSpeedChanges(:);NaNSeparator];
                  jerkDuringField=[jerkDuringField(:); currJerks(:); NaNSeparator];
                  allElapsedTimesInFieldSec=[allElapsedTimesInFieldSec(:);currFieldElapsedTimes(:);NaNSeparator];
                  allElapsedDistsInFieldCm=[allElapsedDistsInFieldCm(:);currFieldElapsedDists(:);NaNSeparator];
                   phaseDiffPerTimeDuringField=[phaseDiffPerTimeDuringField(:);currFieldPhaseDiffs(:);NaNSeparator];
                  allTimesInField=[allTimesInField(:);currFieldTimes(:);NaNSeparator];
                  allCycleDursInField=[allCycleDursInField(:);currCycleDurs(:);NaNSeparator];
                  allCycleNumsInField=[allCycleNumsInField(:);currFieldCycleNums(:);NaNSeparator];
                  allCycleStartTimesInField=[allCycleStartTimesInField(:);currFieldCycleStartTimes(:);NaNSeparator];
                  
                  allPhasesInField=[allPhasesInField(:);currLapPhases(:);NaNSeparator];
                  allLocalSpikeTimesFromFieldEntry=[allLocalSpikeTimesFromFieldEntry(:);currSpikeTimesFromFieldEntry(:);NaNSeparator];
                  
                  allLapNumsInField=[allLapNumsInField(:);currLapNums(:);NaNSeparator];
                  
            end
        
        end
        
        
    end
    
    cycleChangeIdxesInField=[diff(allCycleNumsInField);NaN];
    
    distanceInFieldCmPerTime=distanceFieldFracPerTime;
    meanFieldWidthCm=nanmean(fieldWidthPerLap);
    %meanFieldTimeWidthSec=nanmedian(timeInFieldPerLap);
    %since non-gaussian, use most common field duration (hist peak)
    %meanFieldTimeWidthSec=nanmean(timeInFieldPerLap); 
    [normalTimeLaps,outlierLapNums,outlierLapTimes]=deleteoutliers(timeInFieldPerLap);
    meanFieldTimeWidthSec=nanmean(normalTimeLaps); 
    
    %plotInDelXDelTDelPspace(allElapsedDistsInFieldCm,allElapsedTimesInFieldSec,allPhasesInField,phaseDiffPerTimeDuringField,nanmean(timeInFieldPerLap),allLapNumsInField)
      %[inFieldDelXTP]=plotInDelXDelTDelPspace(allElapsedDistsInFieldCm,allLocalSpikeTimesFromFieldEntry,allPhasesInField,phaseDiffPerTimeDuringField,nanmean(timeInFieldPerLap),nanmean(fieldWidthPerLap),allLapNumsInField);
      %[inFieldDelXTP]=plotInDelCycleXDelCycleTDelCyclePspace(allElapsedDistsInFieldCm,allLocalSpikeTimesFromFieldEntry,allPhasesInField,phaseDiffPerTimeDuringField,nanmean(timeInFieldPerLap),...
      %    nanmean(fieldWidthPerLap),allLapNumsInField,allCycleNumsInField,allTimesInField,absSpeedsDuringField,cycleStartTimes,distanceInFieldCmPerTime,timeAxis);

       [inFieldDelXTP,inFieldDelXTC,inFieldDelXTS,twoCycleTimeDiffPerTime,twoCycleSpaceDiffPerTime]=plotInDelCycleXDelCycleTDelCyclePspace(allElapsedDistsInFieldCm,allLocalSpikeTimesFromFieldEntry,allPhasesInField,phaseDiffPerTimeDuringField,meanFieldTimeWidthSec,...
          nanmean(fieldWidthPerLap),allLapNumsInField,allCycleNumsInField,allTimesInField,absSpeedsDuringField,cycleStartTimes,distanceInFieldCmPerTime,timeAxis,fieldWidthPerLap,timeInFieldPerLap,allCycleDursInField,outlierLapNums);


       isInOutlierLapInField=zeros(size(allLapNumsInField));
      for i=1:length(allLapNumsInField)
          if(min(abs(allLapNumsInField(i)-outlierLapNums))==0)
              isInOutlierLapInField(i)=1;
          end
      end
    
        isInOutlierLapPerTime=false(size(trueLapNumPerTime));
      for i=1:length(trueLapNumPerTime)
          if(min(abs(trueLapNumPerTime(i)-outlierLapNums))==0)
              isInOutlierLapPerTime(i)=true;
          end
      end
      
      
      dispPhasesInField=allPhasesInField;
      dispPhasesInField(logical(isInOutlierLapInField))=NaN;
      dispSpeedsInField=absSpeedsDuringField;
      dispSpeedsInField(logical(isInOutlierLapInField))=NaN;
      dispCycDursInField=allCycleDursInField;
      dispCycDursInField(logical(isInOutlierLapInField))=NaN;
      
      dispPhasesInField(dispSpeedsInField<5)=NaN;
    %{
      figure;
       
        for li=1:length(phasesPerLap)
            plot(phasesPerLap{li},'k.')
            allPhasesInLaps=[allPhasesInLaps;phasesPerLap{li}];
            hold on
        end
        
        dispIdxes=~isnan(allPhasesInLaps) & abs(phaseDiffPerTimeDuringField)>0;
        avgFieldWidth=nanmean(fieldWidthPerLap);
        avgFieldTime=nanmean(timeInFieldPerLap);
        
        avgCycleDur=nanmean(allCycleDursInField);
        %meanSpeed=avgFieldWidth/avgFieldTime;
        %speedToPreserve=avgFieldWidth/avgCycleDur;
        typicalThetaSequenceLength=50; %cm
        meanSpeed=avgFieldWidth/avgFieldTime;
        meanSpeed=avgFieldWidth/1;
        speedToPreserve=typicalThetaSequenceLength/avgCycleDur;
        x=allElapsedDistsInField(dispIdxes);
        t=allElapsedTimesInField(dispIdxes);
    scatter3(x,t,allPhasesInLaps(dispIdxes),30,allPhasesInLaps(dispIdxes),'filled')
        
        [xprime,tprime]=lorentzTransform(x,t,meanSpeed,speedToPreserve);
        figure;
         subplot(2,1,1)
        scatter3(x,t,allPhasesInLaps(dispIdxes),30,allPhasesInLaps(dispIdxes),'filled')
        colormap(gca,jet)
        colorbar
         ylim([0 1])
        zlim([0 360])
        caxis([0 360])
        subplot(2,1,2)
        scatter3(xprime,tprime,allPhasesInLaps(dispIdxes),30,allPhasesInLaps(dispIdxes),'filled')
        colormap(gca,jet)
        colorbar
        %xlim([0 1])
        zlim([0 360])
        caxis([0 360])
        ylim([0 1])
        xlim([-max(xprime) Inf])
        %zlim([-360 0])
    %figure;
    %plot(absSpeedsDuringField, phaseDiffPerTimeDuringField,'k.')
        
        figure;
        histogram(averageStartPhasePerLap,30,'FaceColor','r')
        hold on
        histogram(averageEndPhasePerLap,30,'FaceColor','b')
        
      
        
        estimatedPhaseRange=prctile(averageStartPhasePerLap,95)-prctile(averageEndPhasePerLap,5);
        title(sprintf('%.2f',estimatedPhaseRange))
        %xlim([0 360])
        %}

    
        distanceFieldFracPerTime=distanceFieldFracPerTime/nanmean(fieldWidthPerLap);
    %distanceFieldFracPerTime=distanceFieldFracPerTime/nanmedian(fieldWidthPerLap);
    %timeSinceLastFieldEntry=timeSinceLastFieldEntry/nanmean(timeSinceLastFieldEntry);
     %timeSinceLastFieldEntry=timeSinceLastFieldEntry/nanmedian(timeInFieldPerLap); %has a heavy tail, median better measure of centrality
         timeSinceLastFieldEntry=timeSinceLastFieldEntry/meanFieldTimeWidthSec;

         avgPeakRateAcrossLaps=nanmean(peakSpikeRatePerLap);
         
     
         %{
    figure; 
    subplot(2,2,1)
    histogram(fieldWidthPerLap,15)
    subplot(2,2,2)
    histogram(timeInFieldPerLap,15)
    

    subplot(2,2,3)
    
    plot(timeSinceLastFieldEntryNormLap,distanceFieldFracPerTimeNormLap,'k-')
     subplot(2,2,4)
    
    plot(timeSinceLastFieldEntry,distanceFieldFracPerTime,'r-')
    
    xlim([0 1])
     %}
    

    
    
    %{
    figure; 
  
    %plot(timeAxis,posAlongTrackPerTime,'k-')
    dispTimeAxis=timeAxis;
        dispTimeAxis(~dirIdxes)=NaN;
        
        dispPos=posAlongTrackPerTime;
        dispPos(~dirIdxes)=NaN;
        plot(dispTimeAxis,dispPos,'k-')
    hold on
    for li=1:numLaps
        plot(fieldEntryTime(li),fieldEntryPos(li),'bo')
        plot(fieldExitTime(li),fieldExitPos(li),'ro')
    end
    plot(xlim,[manualFieldStartCm manualFieldStartCm],'b-')
    plot(xlim,[manualFieldEndCm manualFieldEndCm],'r-')
  
    
    if(manualSelectedDI==2)
        disp('')
        
    end
   %}
    %only consider speed change per position change during earnest laps
  
  
   
    
    speedChangePerPosChangePerLap=NaN(numLaps,maxTimeIdxDuringLap);
    
    
    for li=1:numLaps
        
        currLapTimeIdxes=trueLapNumPerTime==li;
        currLapSpeedChanges=speedChangePerPosChangePerTime(currLapTimeIdxes);
        %reverse direction travel
        if(providedLapStartPosition(li)>trackLengthCm/2)
            currLapSpeedChanges=-currLapSpeedChanges;
        end
        speedChangePerPosChangePerLap(li,1:length(currLapSpeedChanges))=currLapSpeedChanges;
    end
    
   
    [~,durationSortedLapIdxes]=sort(trueLapDurations);
    durationSortedLapSpeedChanges=speedChangePerPosChangePerLap(durationSortedLapIdxes,:);
    %figure; imagesc(durationSortedLapSpeedChanges)
    %figure;plot(durationSortedLapSpeedChanges(1:100,:)')
    
    spaceTimeNumBins=30;
     %spaceTimeNumBins=25;
    %maxLapTime=7;
      maxLapTime=5;
      %maxLapTime=3.5;
      
      %minFieldTime=0.5;
    timeSinceLastTrueLapStartPerTime(timeSinceLastTrueLapStartPerTime>maxLapTime)=NaN;
    timeSinceLastFieldEntry(timeSinceLastFieldEntry>maxFieldTime)=NaN;
    %timeSinceLastFieldEntry(timeSinceLastFieldEntry<minFieldTime)=NaN;
    %stH=figure 
    %for di=1:2
        youShouldFlipManualPlaceBounds=0;
        
        if(manualSelectedDI==1)
            distanceFromStartPerTime=trackLengthCm-posAlongTrackPerTime;
            
        else
            distanceFromStartPerTime=posAlongTrackPerTime;
            %youShouldFlipManualPlaceBounds=1;
        end
        [distBins,timeBins,meanRateSmooth,dataCountPerBin] = makeHeatMapAvgOf3rdVar(distanceFromStartPerTime(dirIdxes),timeSinceLastTrueLapStartPerTime(dirIdxes),localSpikeRatePerTime(dirIdxes),spaceTimeNumBins)
            %subplot(2,1,di)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %____Take out outlier laps - Dec 22,2020
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            dirIdxes=(~isInOutlierLapPerTime) & dirIdxes(:) & ~thisLapHasBackwardsMovementPerTime(:);
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %^^^^Take out outlier laps - Dec 22,2020
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
         placeCellDistPerTime=distanceFromStartPerTime(dirIdxes);
         placeCellTimeInLapPerTime=timeSinceLastTrueLapStartPerTime(dirIdxes);
         
         placeCellTimeInFieldPerTime=timeSinceLastFieldEntry(dirIdxes);
         placeCellFieldFracPerTime=distanceFieldFracPerTime(dirIdxes);
         
         placeCellRatePerTime=localSpikeRatePerTime(dirIdxes);
         placeCellPhasePerTime=localSpikePhasePerTime(dirIdxes);
         placeCellThetaAmpPerTime=localCycleAmpPerTime(dirIdxes);
         placeCellThetaDurPerTime=localCycleDurPerTime(dirIdxes);
         
         %no spikes means nan phase
         placeCellPhasePerTime(placeCellPhasePerTime==0)=NaN;
         
         placeCellSpeedPerTime=absSpeedAlongTrackPerTime(dirIdxes);
         
         %[peakRateSpaceTime,peakRateSpaceTimeIdx]=max(meanRateSmooth(:));        
         %[peakDistBin, peakTimeBin]=ind2sub(size(meanRateSmooth),peakRateSpaceTimeIdx);
         %peakDist=distBins(peakDistBin);
         %peakTime=timeBins(peakTimeBin);
         
         goodRegionDist=distBins>manualFieldStartCm & distBins<manualFieldEndCm;
         goodRegionRate=meanRateSmooth;
         goodRegionRate(~goodRegionDist,:)=NaN;
         
         manualMaxTime=3;
         goodRegionTime=timeBins<manualMaxTime;
         goodRegionRate(:,~goodRegionTime)=NaN;
         
         spacetimerateCOM=centerOfMass(goodRegionRate);
         try
         comDist=distBins(round(spacetimerateCOM(1)));
         comTime=timeBins(round(spacetimerateCOM(2)));
         catch
             comDist=0;
             comTime=0;
         end
         if(zeroFieldToCOM)
            placeCellDistPerTime=placeCellDistPerTime-comDist;
             placeCellTimeInLapPerTime=placeCellTimeInLapPerTime-comTime;
         else
             comDist=0;
             comTime=0;
         end
         %{
        omarPcolor(distBins-comDist,timeBins-comTime,meanRateSmooth',stH)
        xlabel('Distance from field COM (cm)')
        ylabel('Time in true lap from field COM (sec)')
        colormap(jet)
        colorbar
        %title(sprintf('Direction %d',di))
        disp('')
         %}
     %end
    %makeHeatMapAvgOf3rdVar(posAlongTrackPerTime(dirIdxes),absSpeedAlongTrackPerTime(dirIdxes),localSpikeRatePerTime(dirIdxes),spaceTimeNumBins)
        %makeHeatMapAvgOf3rdVar(posAlongTrackPerTime(dirIdxes),speedChangePerPosChangePerTime(dirIdxes),localSpikeRatePerTime(dirIdxes),spaceTimeNumBins)

        
            
        nonNanPhaseIDs=~isnan(placeCellPhasePerTime);
        
            %scaledTime=scaledata(0,1,placeCellTimeInLapPerTime);
            %scaledDist=scaledata(0,1,placeCellDistPerTime);
            
            findBestStartShiftPhases=placeCellPhasePerTime(nonNanPhaseIDs);
            
            %scaledTime=placeCellTimeInLapPerTime(nonNanPhaseIDs)+comTime;
            scaledTime=placeCellTimeInFieldPerTime(nonNanPhaseIDs)+comTime;
            
            scaledDist=placeCellDistPerTime(nonNanPhaseIDs)+comDist;
            scaledTime=scaledTime(:);
            scaledDist=scaledDist(:);
            
           
        
            
            avgValuePerShift=NaN(360,1);
            
            scaledTime=zscoreLFP(scaledTime);
            scaledDist=zscoreLFP(scaledDist);
            
            scaledTime=scaledTime-nanmin(scaledTime);
            scaledDist=scaledDist-nanmin(scaledDist);
             %boundaryScaledTime=0.4*range(scaledTime);
            %boundaryScaledDist=0.4*range(scaledDist);
            
             %boundaryScaledTime=0.5*range(scaledTime);
            %boundaryScaledDist=0.5*range(scaledDist);
                %startingRegion=scaledTime<boundaryScaledTime & scaledDist<boundaryScaledDist;
                
                spaceTimeOriginDist=sqrt(scaledTime.^2+scaledDist.^2);
                %boundarySTdist=prctile(spaceTimeOriginDist,30);
                boundarySTdist=prctile(spaceTimeOriginDist,20);
                boundarySTdistEnd=prctile(spaceTimeOriginDist,80); %bigger area farther out
                 %boundarySTdist=prctile(spaceTimeOriginDist,10);
                 
                 %boundarySTdist=prctile(spaceTimeOriginDist,20);
                 
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 %standardize circ shift according to spacetime distnace start
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 %boundarySTdist=2; 
                 %boundarySTdist=prctile(spaceTimeOriginDist,5)
                 
                 %figure; histogram(spaceTimeOriginDist)
                startingRegion=spaceTimeOriginDist<=boundarySTdist;
                
                endingRegion=spaceTimeOriginDist>=boundarySTdistEnd;
            
            bestShiftDeg=0;
            for si=1:length(avgValuePerShift)
                startingRegionPhases=findBestStartShiftPhases(startingRegion)+si;
                startingRegionPhases=mod(startingRegionPhases,360);
                
                endingRegionPhases=findBestStartShiftPhases(endingRegion)+si;
                endingRegionPhases=mod(endingRegionPhases,360);
                
               avgValuePerShift(si)=nanmean(startingRegionPhases)-nanmean(endingRegionPhases); %less penalty farther out?
            end
            [~,bestShiftDeg]=nanmax(avgValuePerShift);
            %}
          
            %bestShiftDeg=data.manualPrecessionStartPhaseDeg-360;
            
            %[bestShiftDeg]=getBestLinearCorrShift(findBestStartShiftPhases,sqrt(scaledTime.^2+scaledDist.^2)); %circ shift to corr with distance from origin in space time
          
          placeCellPhasePerTime=mod(placeCellPhasePerTime+bestShiftDeg,360);
          
          
          
          %placeCellPhasePerTime= mod(placeCellPhasePerTime-(manualPrecessionStartPhaseDeg-360),360);
          
          %placeCellPhasePerTime=scaledata(placeCellPhasePerTime,0, 360);
        
    save(filePath,'trueLapStartTimes','trueLapEndTimes','trueLapDurations', 'trueLapStartPositions','trueLapEndPositions','durationSortedLapIdxes', 'speedChangePerPosChangePerLap',...
        'posAlongTrackPerTime','absSpeedAlongTrackPerTime','placeCellSpeedPerTime','trueLapNumPerTime','timeSinceLastTrueLapStartPerTime','placeCellRatePerTime','placeCellPhasePerTime', ...
        'placeCellDistPerTime','placeCellTimeInLapPerTime','maxLapTime','placeCellThetaDurPerTime','placeCellThetaAmpPerTime','comDist','comTime','fieldEntryTime','fieldExitTime','placeCellTimeInFieldPerTime',...
        'speedChangePerPosChangePerTime','youShouldFlipManualPlaceBounds','trackLengthCm','manualSelectedDI','timeSinceLastFieldEntry','placeCellFieldFracPerTime','avgPeakRateAcrossLaps','peakSpikeRatePerLap', ...
        'avgWaveform','inFieldDelXTP','allElapsedDistsInFieldCm','dispSpeedsInField','allElapsedTimesInFieldSec','meanFieldTimeWidthSec','inFieldDelXTC','inFieldDelXTS','dispPhasesInField','cycleChangeIdxesInField',...
        'dispCycDursInField','meanFieldWidthCm','twoCycleTimeDiffPerTime','twoCycleSpaceDiffPerTime','speedChangesPerPosDuringField','jerkDuringField','-append')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %phase portrait movie
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    showMovie=0;
    %providedLapStartTimes
    if(showMovie)
        figure
        for stopDisp=1:12000
           
            currTime=timeAxis(stopDisp);
      
            [timeFromLapStartClark, currLapNumClark]=min(abs(currTime-providedLapStartTimes));
            [timeFromLapStart, currLapNum]=min(abs(currTime-trueLapStartTimes));
            [timeFromLapEnd, currLapNumEnd]=min(abs(currTime-trueLapEndTimes));
            
            dispNewLapNum=0;
            if(timeFromLapStart<0.03)
                dispNewLapNum=1;
            end
             dispNewLapNumClark=0;
            if(timeFromLapStartClark<0.03)
                dispNewLapNumClark=1;
            end
            
            dispNewLapNumEnd=0;
            if(timeFromLapEnd<0.03)
                dispNewLapNumEnd=1;
            end
        subplot(3,1,1)
        %if(dispNewLapNum==1 && dispNewLapNumEnd==0)
        %p1=plot(timeAxis(1:stopDisp),speedChangePerPosChangePerTime(1:stopDisp),'k-');
        %hold on
        %p2=plot(timeAxis(1:stopDisp),speedChangePerPosChangePerTime(1:stopDisp),'ko','MarkerSize',2);
        currLapID=trueLapNumPerTime(stopDisp);

        if(~isnan(currLapID))
            nonNaNSpeedChange=speedChangePerPosChangePerLap(currLapID,:);
            nonNaNSpeedChange(isnan(nonNaNSpeedChange))=[];
            currLapTimeAxis=(1:length(nonNaNSpeedChange))*approxTimeStep;
            p1=plot(currLapTimeAxis,nonNaNSpeedChange,'k-');
            hold on
            p2=plot(currLapTimeAxis,nonNaNSpeedChange,'ko','MarkerSize',2);
        end

        %end
        %p1=plot(posAlongTrackPerTime(1:stopDisp),(speedAlongTrackPerTime(1:stopDisp)),'k-');
        %hold on
         %p2=plot(posAlongTrackPerTime(1:stopDisp),(speedAlongTrackPerTime(1:stopDisp)),'ko','MarkerSize',2);
          %pDot1=plot(posAlongTrackPerTime(stopDisp),speedAlongTrackPerTime(stopDisp),'bo','MarkerSize',5);
        %plot(posAlongTrackPerTime,speedAlongTrackPerTime,'bo')
        %xlabel('Distance along track (cm)')
        %ylabel('Velocity (cm/s)')
        ylabel('Speed change per cm (cm/s per cm)')
        xlabel('Time (sec)')
        
        subplot(3,1,2)
        p3=plot(posXPerTimeStep(1:stopDisp),posYPerTimeStep(1:stopDisp),'k-');
        hold on
        pDot2=plot(posXPerTimeStep(stopDisp),posYPerTimeStep(stopDisp),'bo','MarkerSize',5);
        xlabel('X distance (cm)')
        ylabel('Y distance (cm)')
        
        subplot(3,1,3)
        p4=plot(timeAxis(1:stopDisp),posAlongTrackPerTime(1:stopDisp),'k-');
        hold on
        
        xlabel('Time (sec)')
        ylabel('Distance along track (cm)')
        if(dispNewLapNum)
            title(sprintf('Start lap %d',currLapNum))
            p5=plot([currTime currTime],[0 120],'k-','LineWidth',3);
            hold on
            %pause(0.5)
        end
        if(dispNewLapNumEnd)
            title(sprintf('End lap %d',currLapNumEnd))
            p5=plot([currTime currTime],[0 120],'r-','LineWidth',3);
            hold on
            %pause(0.5)
        end
        if(dispNewLapNumClark)
            %title(sprintf('Clark start lap %d',currLapNumClark))
            %p5=plot([currTime currTime],[0 120],'Color',getGrayRGB(),'LineWidth',1);
            %hold on
            %pause(0.5)
        end
        
        
        drawnow
        
        %delete(pDot1)
        try
            delete(pDot2)
        end
        try
            delete(p1)
        end
        try
            delete(p2)
        end
        try
            delete(p3)
        end
        try
            delete(p4)
        end
        
        %delete(p5)

        end
    end
    %plot(posAlongTrackPerTime(1:stopDisp),speedPerTimeStep(1:stopDisp),'r-')
    
    %figure
    %plot(posAlongTrackPerTime,(speedAlongTrackPerTime(:)-speedPerTimeStep(:)),'k.')
    
    %for di=1:(length(posAlongTrackPerTime)-1)
        
    %end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Vector projected position is better measure of position along track
    %than either one alone
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %{
    figure;
    subplot(3,1,1)
    plot(posXPerTimeStep,posYPerTimeStep,'k-')
    hold on
    %plot(posXPerTimeStep,posYPerTimeStep,'ko')
    subplot(3,1,2)
    figure
    plot(timeAxis,posAlongTrackPerTime,'k-')
    hold on
     plot(timeAxis,posXPerTimeStep,'r-')
     plot(timeAxis,posYPerTimeStep,'b-')
    %}
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Vector projected speed is better measure of speed along track
    %than provided speed!
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %{
    figure
    %yyaxis left
    normVecSpeed=scaledata(abs(speedAlongTrackPerTime),0,1);
    normXSpeed=scaledata(abs(velocityXYPerTimeStep),0,1);
    
    plot(timeAxis,normVecSpeed,'k-')
    speedPerTimeStep=originalData.frames(withinTaskFrameIdxes,5);
     noiseSpeedThresh=200;
    speedPerTimeStep(speedPerTimeStep>noiseSpeedThresh)=NaN;
    normProvidedSpeed=scaledata(speedPerTimeStep,0,1);
    %yyaxis right
    hold on
    plot(timeAxis,normProvidedSpeed,'r-')
     %plot(timeAxis,normXSpeed,'b-')
    %}
    
    disp('')
    
    
    

function [spikePerCycleInfo] = getPerCycleSpikeInfoClarkFile(fileName,cellNum)
close all
tic
recompute=0;


%fileName='LEM3216_S20190726184722.mat';
%fileName='LEM3246_S20190715190235.mat';
fileRootName=fileName(1:(end-4));
dataDir='/Users/tibinjohn/Downloads/processed_data';
matchingFilePaths=getFilePathsRegex(dataDir,[fileRootName '*.mat']);

%data=load(fullfile(dataDir,fileName));
data=load(matchingFilePaths{1});
%cellNum=3;
%cellNum=11;

taskStartTime=data.events(1,1);
taskEndTime=data.events(2,1); %2nd row of 1st col (linear track end)

speedOverWholeSession=data.binned_vel{1}; %cm/sec, 200msec bins
binWidthSpeed=0.2; %200msec
%speedTimeAxis=linspace(taskStartTime,taskEndTime,length(speedOverWholeSession)+1);
speedTimeAxis=(taskStartTime+binWidthSpeed/2):binWidthSpeed:(taskEndTime-binWidthSpeed/2);
stopLater=0;
if(length(speedTimeAxis)>length(speedOverWholeSession))
    speedTimeAxis=speedTimeAxis(1:length(speedOverWholeSession)); %assume start is correct
    stopLater=1;
end


trackLengthCm=data.maze_size_cm(1);
trackLengthBins=trackLengthCm/3; %3 cm bins

cellIDstr=sprintf('%s_Cell%d',fileRootName,cellNum);

spikeTimes=data.Spikes{cellNum};

%rateTimeBinWidth=0.1; %vs more than theta cycle?
rateTimeBinWidth=0.3; %vs more than theta cycle?
rateTimeBinWidth=0.5; %vs more than theta cycle?
rateTimeBinWidth=0.05; %vs more than theta cycle?
[rateTimeBinCenters,firingRateOverTime]=spikeTimesToFiringRate(spikeTimes,rateTimeBinWidth);

correspondingTetStr=data.spikesID.TetrodeNum{cellNum};
tetNum=str2num(getSubstrBtwnSubstrs(correspondingTetStr,'TT','.mat'));

channelToTetNum=data.lfp.channel_list.tetrode_num;
chNums=find(channelToTetNum==tetNum);

avgWaveformsPerCh=data.avgwave{cellNum};
[~,maxChInd]=max(max(avgWaveformsPerCh'));

offset=(maxChInd-1);
%not tetrode condition
if(max(data.lfp.channel_list.tetrode_num)==length(data.lfp.channel_list.tetrode_num))
    offset=0;
end
chWithGreatestAmp=min(chNums)+offset;

peakFieldMaxRate=-1;

timeAxis=data.lfp.ts(:);
linearTrackIDsOnly=(timeAxis<taskEndTime);

correspondingLFP=data.lfp.signal(chWithGreatestAmp,linearTrackIDsOnly);
timeAxis=timeAxis(linearTrackIDsOnly);

maxTimeLFP=max(timeAxis);
Fs=1/median(diff(timeAxis));

%lowFreq=3;
lowFreq=4;
highFreq=12;
lowFreq=6;
highFreq=12;
%[cycleInfo,cycleTimeAxis,cycleLFP] = rawLFPtoCycleInfo(timeAxis,correspondingLFP,lowFreq,highFreq);
%figure; plot(cycleTimeAxis,cycleLFP)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%UPSAMPLE LFP TO GET MORE PRECISE CYCLE TIMES - PREVIOUSLY 5ms precision now 0.5msec precision
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
interpFact=10;
interpIdxs=1:(1/interpFact):length(timeAxis);
interpolatedTimeAxis=interp1(1:length(timeAxis),timeAxis,interpIdxs,'linear');

interpLFP=interp1(timeAxis,correspondingLFP,interpolatedTimeAxis,'linear');

%[cycleInfo] = ...
%        getCycleInfoLFP(timeAxis, correspondingLFP, Fs, lowFreq, highFreq, cellIDstr,0);
[cycleInfo] = ...
        getCycleInfoLFP(interpolatedTimeAxis,interpLFP, Fs*interpFact, lowFreq, highFreq, cellIDstr,0);

minTimes=cycleInfo.cycleMinTimes;
maxTimes=cycleInfo.cycleMaxTimes;

%thetaAmpZPerCycle=zscoreLFP(cycleInfo.cycleMaxAmp-cycleInfo.cycleMinAmp);
thetaAmpZPerCycle=(cycleInfo.cycleMaxAmp-cycleInfo.cycleMinAmp);

numCycles=length(cycleInfo.cycleMaxTimes);
[spikeAssignMax spikeAssignMin spikeAssignZero] = assignSpikesToCycles2017(spikeTimes, minTimes, maxTimes);

mode=2; %phase from max to max
mode=1;%phase from startMin to endMin
calcTimeOffsets=1;
[spikePhase spikeOffsets spikeOffsetsEnd] = assignPhaseToSpikes(spikeTimes, spikeAssignMax, spikeAssignMin, minTimes, maxTimes, mode, calcTimeOffsets, spikeAssignZero);

[spikePerCycleInfo]=getFirstSpikePhasePerCycle(cycleInfo,spikePhase,spikeAssignMax);

lapStartTimes=[data.linear_track{1}.lapinfo.laps.start_ts];
lapDirections=[data.linear_track{1}.lapinfo.laps.direction];

%data.linear_track{1}.right{3}.fields{1}
%fieldStartPos=
%fieldEndPos=

posTimeAxis=data.linear_track{1}.nonlinearFrames(:,1);
posPerTimeStep=data.linear_track{1}.nonlinearFrames(:,2);
posPerTimeStepNorm=scaledata(posPerTimeStep,0,1);

posTimeAxisStep=median(diff(posTimeAxis));

timeSpentPosBinSize=0.01;
timeSpentBinEdges=0:timeSpentPosBinSize:1; %norm distance
timeSpentAtEachPosition=histcounts(posPerTimeStepNorm,timeSpentBinEdges);
timeSpentAtEachPosition=timeSpentAtEachPosition*posTimeAxisStep;

posPerCycle=NaN(numCycles,1);
distTravelledPerCycleNormTrack=NaN(numCycles,1);
normPosPerCycle=NaN(numCycles,1);
elapsedTimeInLapPerCycle=NaN(numCycles,1);
lapDirectionPerCycle=NaN(numCycles,1);
lapNumPerCycle=NaN(numCycles,1);
isInFieldPerCycle=zeros(numCycles,1);

speedPerCycle=NaN(numCycles,1);

peakFieldStartPosNorm=NaN(2,1);
peakFieldEndPosNorm=NaN(2,1);

firingRatePerCycle=NaN(numCycles,1);

for i=1:numCycles
    cycleStartTime=maxTimes(i);
    if(i<numCycles)
        cycleEndTime=maxTimes(i+1);
    else
        cycleEndTime=maxTimeLFP;
    end
    posCycleIdxes=find(cycleStartTime<=posTimeAxis & posTimeAxis<=cycleEndTime);
    rateCycleIdxes=find(cycleStartTime<=rateTimeBinCenters & rateTimeBinCenters<=cycleEndTime);
    if(~isempty(rateCycleIdxes))
	firingRatePerCycle(i)=nanmean(firingRateOverTime(rateCycleIdxes));
    end

    %posPerCycle(i)=mean(posPerTimeStep(posCycleIdxes)); %distance from within cycle variation added in later by phase proportion
    %currPosNorm=mean(posPerTimeStepNorm(posCycleIdxes));
    
    if(~isempty(posCycleIdxes)) 
        posPerCycle(i)=min(posPerTimeStep(posCycleIdxes)); %distance from within cycle variation added in later by phase proportion
        currPosNorm=min(posPerTimeStepNorm(posCycleIdxes));
        normPosPerCycle(i)=currPosNorm;
   
        distTravelledPerCycleNormTrack(i)=range(posPerTimeStepNorm(posCycleIdxes));
    end
    speedCycleIdxes=find(cycleStartTime<=speedTimeAxis & speedTimeAxis<cycleEndTime+binWidthSpeed);%guarantee at least one value
    speedPerCycle(i)=nanmean(speedOverWholeSession(speedCycleIdxes)); 
    
    [~,latestLapStartIdx]=min(abs(cycleStartTime-lapStartTimes));
    latestLapStartTime=lapStartTimes(latestLapStartIdx);
    if(latestLapStartTime>cycleStartTime)
        if(latestLapStartIdx>1)
            latestLapStartIdx=latestLapStartIdx-1;
            %latestLapStartTime=lapStartTimes(latestLapStartIdx-1);
            latestLapStartTime=lapStartTimes(latestLapStartIdx);
        end
    end
    
    elapsedTimeInLapPerCycle(i)=cycleStartTime-latestLapStartTime; %time from within cycle variation added in later by phase proportion

    lapDirectionPerCycle(i)=lapDirections(latestLapStartIdx);
    lapNumPerCycle(i)=latestLapStartIdx;
    if(lapDirectionPerCycle(i)<0 || isnan(lapDirectionPerCycle(i)))
        thisCellFields=data.linear_track{1}.left{cellNum}.fields;
        dirID=1;
    else
        thisCellFields=data.linear_track{1}.right{cellNum}.fields;
        dirID=2;
    end
    
   thisDirMaxFieldRates=NaN(length(thisCellFields),1);
    %choose only max firing rate field for this direction
    for fi=1:length(thisCellFields)
        thisDirMaxFieldRates(fi)=thisCellFields{fi}.peakFR;
    end
    
    [~,maxFi]=max(thisDirMaxFieldRates);
    %for fi=1:length(thisCellFields)
        currFieldStartNorm=thisCellFields{maxFi}.start/trackLengthBins;
        currFieldEndNorm=thisCellFields{maxFi}.stop/trackLengthBins;
        
        maxFiringRate=thisCellFields{maxFi}.peakFR;
        if(maxFiringRate>peakFieldMaxRate)
             peakFieldStartPosNorm(dirID)=currFieldStartNorm;
             peakFieldEndPosNorm(dirID)=currFieldEndNorm;
             peakFieldMaxRate(dirID)=maxFiringRate;
        end
        
        if(~isempty(posCycleIdxes) && currPosNorm>currFieldStartNorm && currPosNorm<currFieldEndNorm)
            isInFieldPerCycle(i,dirID)=1;
        end
    %end
    
    
end %loop over theta cycles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%store variable values per theta cycle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%isInFieldPerCycle(isInFieldPerCycle

spikePerCycleInfo.filePath=fullfile(dataDir,fileName);
spikePerCycleInfo.cellNum=cellNum;
spikePerCycleInfo.LFPlowFreq=lowFreq;
spikePerCycleInfo.LFPhighFreq=highFreq;

spikePerCycleInfo.lapStartTimes=lapStartTimes;

spikePerCycleInfo.lapDirections=lapDirections;
spikePerCycleInfo.allSpikeTimes=spikeTimes;
spikePerCycleInfo.posPerCycle=posPerCycle;
spikePerCycleInfo.distTravelledPerCycleNormTrack=distTravelledPerCycleNormTrack;

spikePerCycleInfo.thetaAmpZPerCycle=thetaAmpZPerCycle;
spikePerCycleInfo.speedPerCycle=speedPerCycle;
spikePerCycleInfo.durationPerCycleSec=cycleInfo.cycleDurations;

spikePerCycleInfo.firingRatePerCycle=firingRatePerCycle;
spikePerCycleInfo.firingRateOverTime=firingRateOverTime;
spikePerCycleInfo.rateTimeBinCenters=rateTimeBinCenters;

spikePerCycleInfo.elapsedTimeInLapPerCycle=elapsedTimeInLapPerCycle;
spikePerCycleInfo.lapDirectionPerCycle=lapDirectionPerCycle;
spikePerCycleInfo.lapNumPerCycle=lapNumPerCycle;
spikePerCycleInfo.isInFieldPerCycle=isInFieldPerCycle;
spikePerCycleInfo.cellIDstr=cellIDstr;
spikePerCycleInfo.normPosPerCycle=normPosPerCycle;
spikePerCycleInfo.posPerCycleCm=normPosPerCycle*trackLengthCm;

spikePerCycleInfo.peakFieldStartPosNorm=peakFieldStartPosNorm;
spikePerCycleInfo.peakFieldEndPosNorm=peakFieldEndPosNorm;
spikePerCycleInfo.peakFieldMaxRate=peakFieldMaxRate;
spikePerCycleInfo.meanSpeed=data.BasicLoco.MeanVelocity(1); %linear track mean speed
spikePerCycleInfo.taskEndTime=taskEndTime;
spikePerCycleInfo.timeSpentAtEachPosition=timeSpentAtEachPosition;

if(stopLater)
    %figure;  plot(zscoreLFP(spikePerCycleInfo.posPerCycle)*20);hold on; plot(spikePerCycleInfo.speedPerCycle);
    %fds
end
%see data.linear_track{1}.right{3}.fields{1}.ThPrecess
%posPerCycle=

toc

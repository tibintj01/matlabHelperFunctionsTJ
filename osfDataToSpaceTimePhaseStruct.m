function [spikePerCycleInfo] = getPerCycleSpikeInfoClarkFile(fileName,cellNum)
close all; clear all; clc
tic

%fileName='LEM3216_S20190726184722.mat';
fileName='LEM3246_S20190715190235.mat';
fileRootName=fileName(1:(end-4));
dataDir='/Users/tibinjohn/Downloads/processed_data';
data=load(fullfile(dataDir,fileName));
cellNum=3;
cellNum=11;

taskEndTime=data.events(2,1); %2nd row of 1st col (linear track end)


trackLengthCm=data.maze_size_cm(1);
trackLengthBins=trackLengthCm/3; %3 cm bins

cellIDstr=sprintf('%s_Cell%d',fileRootName,cellNum);

spikeTimes=data.Spikes{cellNum};

correspondingTetStr=data.spikesID.TetrodeNum{cellNum};
tetNum=str2num(getSubstrBtwnSubstrs(correspondingTetStr,'TT','.mat'));

channelToTetNum=data.lfp.channel_list.tetrode_num;
chNums=find(channelToTetNum==tetNum);

avgWaveformsPerCh=data.avgwave{cellNum};
[~,maxChInd]=max(max(avgWaveformsPerCh'));

chWithGreatestAmp=min(chNums)+(maxChInd-1);
peakFieldMaxRate=-1;

timeAxis=data.lfp.ts(:);
linearTrackIDsOnly=(timeAxis<taskEndTime);

correspondingLFP=data.lfp.signal(chWithGreatestAmp,linearTrackIDsOnly);
timeAxis=timeAxis(linearTrackIDsOnly);

maxTimeLFP=max(timeAxis);
Fs=1/median(diff(timeAxis));

lowFreq=3;
highFreq=12;
%[cycleInfo,cycleTimeAxis,cycleLFP] = rawLFPtoCycleInfo(timeAxis,correspondingLFP,lowFreq,highFreq);
%figure; plot(cycleTimeAxis,cycleLFP)
[cycleInfo] = ...
        getCycleInfoLFP(timeAxis, correspondingLFP, Fs, lowFreq, highFreq, cellIDstr,0);

minTimes=cycleInfo.cycleMinTimes;
maxTimes=cycleInfo.cycleMaxTimes;

numCycles=length(cycleInfo.cycleMaxTimes);
[spikeAssignMax spikeAssignMin spikeAssignZero] = assignSpikesToCycles2017(spikeTimes, minTimes, maxTimes);

mode=2; %phase from max to max
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



posPerCycle=NaN(numCycles,1);
normPosPerCycle=NaN(numCycles,1);
elapsedTimeInLapPerCycle=NaN(numCycles,1);
lapDirectionPerCycle=NaN(numCycles,1);
lapNumPerCycle=NaN(numCycles,1);
isInFieldPerCycle=zeros(numCycles,1);

peakFieldStartPosNorm=NaN(2,1);
peakFieldEndPosNorm=NaN(2,1);

for i=1:numCycles
    cycleStartTime=maxTimes(i);
    if(i<numCycles)
        cycleEndTime=maxTimes(i+1);
    else
        cycleEndTime=maxTimeLFP;
    end
    posCycleIdxes=find(cycleStartTime<posTimeAxis & posTimeAxis<cycleEndTime);
    posPerCycle(i)=mean(posPerTimeStep(posCycleIdxes));
    currPosNorm=mean(posPerTimeStepNorm(posCycleIdxes));
    normPosPerCycle(i)=currPosNorm;
    
    [~,latestLapStartIdx]=min(abs(cycleStartTime-lapStartTimes));
    latestLapStartTime=lapStartTimes(latestLapStartIdx);
    if(latestLapStartTime>cycleStartTime)
        if(latestLapStartIdx>1)
            latestLapStartIdx=latestLapStartIdx-1;
            %latestLapStartTime=lapStartTimes(latestLapStartIdx-1);
            latestLapStartTime=lapStartTimes(latestLapStartIdx);
        end
    end
    
    elapsedTimeInLapPerCycle(i)=cycleStartTime-latestLapStartTime;
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
        
        if(currPosNorm>currFieldStartNorm && currPosNorm<currFieldEndNorm)
            isInFieldPerCycle(i,dirID)=1;
        end
    %end
    
    
end

%isInFieldPerCycle(isInFieldPerCycle

spikePerCycleInfo.filePath=fullfile(dataDir,fileName);
spikePerCycleInfo.cellNum=cellNum;
spikePerCycleInfo.LFPlowFreq=lowFreq;
spikePerCycleInfo.LFPhighFreq=highFreq;

spikePerCycleInfo.lapStartTimes=lapStartTimes;

spikePerCycleInfo.lapDirections=lapDirections;
spikePerCycleInfo.allSpikeTimes=spikeTimes;
spikePerCycleInfo.posPerCycle=posPerCycle;
spikePerCycleInfo.elapsedTimeInLapPerCycle=elapsedTimeInLapPerCycle;
spikePerCycleInfo.lapDirectionPerCycle=lapDirectionPerCycle;
spikePerCycleInfo.lapNumPerCycle=lapNumPerCycle;
spikePerCycleInfo.isInFieldPerCycle=isInFieldPerCycle;
spikePerCycleInfo.cellIDstr=cellIDstr;
spikePerCycleInfo.normPosPerCycle=normPosPerCycle;
spikePerCycleInfo.peakFieldStartPosNorm=peakFieldStartPosNorm;
spikePerCycleInfo.peakFieldEndPosNorm=peakFieldEndPosNorm;
spikePerCycleInfo.peakFieldMaxRate=peakFieldMaxRate;
spikePerCycleInfo.meanSpeed=data.BasicLoco.MeanVelocity(1); %linear track mean speed
spikePerCycleInfo.taskEndTime=taskEndTime;


%see data.linear_track{1}.right{3}.fields{1}.ThPrecess
%posPerCycle=

toc
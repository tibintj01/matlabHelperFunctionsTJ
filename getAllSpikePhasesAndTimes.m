function [allSpikeTimesWithPhases]=getAllSpikePhasesAndTimes(originalData,data)

    cellNum=data.spikePerCycleInfo.cellNum;
    spikeTimes=originalData.Spikes{cellNum};
    
    
    fileRootName=getFileNameFromPath(data.spikePerCycleInfo.filePath);
    cellIDstr=sprintf('%s_Cell%d',fileRootName,cellNum);
    
    taskStartTime=originalData.events(1,1);
    taskEndTime=originalData.events(2,1); %2nd row of 1st col (linear track end)

    correspondingTetStr=originalData.spikesID.TetrodeNum{cellNum};
    tetNum=str2num(getSubstrBtwnSubstrs(correspondingTetStr,'TT','.mat'));

    channelToTetNum=originalData.lfp.channel_list.tetrode_num;
    chNums=find(channelToTetNum==tetNum);

    avgWaveformsPerCh=originalData.avgwave{cellNum};
    [~,maxChInd]=max(max(avgWaveformsPerCh'));

    offset=(maxChInd-1);
    %not tetrode condition
    if(max(originalData.lfp.channel_list.tetrode_num)==length(originalData.lfp.channel_list.tetrode_num))
        offset=0;
    end
    chWithGreatestAmp=min(chNums)+offset;


    timeAxis=originalData.lfp.ts(:);
    linearTrackIDsOnly=(timeAxis<=taskEndTime);

    correspondingLFP=originalData.lfp.signal(chWithGreatestAmp,linearTrackIDsOnly);
    timeAxis=timeAxis(linearTrackIDsOnly);

    maxTimeLFP=max(timeAxis);
    Fs=1/median(diff(timeAxis));

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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %get cycle information
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %save filtered LFP
    saveFilteredLFP=0;
    [cycleInfo] = ...
            getCycleInfoLFP(interpolatedTimeAxis,interpLFP, Fs*interpFact, lowFreq, highFreq, cellIDstr,saveFilteredLFP);
    


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
    
    allSpikeTimesWithPhases=[spikeTimes(:), spikePhase(:)];

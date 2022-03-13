%% HELPER FILE FOR human_spike_LFP_relationship_central
%  Just moved some variable initiation (setting to NaN) code here
%  Make things cleaner by saving data in spikeNS5Struct, cycleStruct, spikePhaseStruct, spikeRateCycleStruct, spikeRateSpecWindowStruct

%% spikeNS5Struct
spikeNS5Struct.spikeNS5FileInd = NaN(1, numSpikes);
spikeNS5Struct.spikeTimesWithinNS5File = NaN(1,numSpikes); % adjusted to start from 0 at the begnning of each 40 min NS5 file
spikeNS5Struct.ns5SpikeMap = NaN(dataInfo(currSes).numNS5Files, 2);

%%  cycleStruct: use for cycle->rate and cycle->amp/duration/slope and cycle->spectralPower information
%   cycleStruct is stored separately - one file per LFP (since it takes up
%   60 MB storage per cell/LFP combintation otherwise - just load it when
%   loading the cell/LFP information
cycleStruct.numCyclesPerNS5File = NaN(1,dataInfo(currSes).numNS5Files);
cycleStruct.ns5CycleMap = NaN(dataInfo(currSes).numNS5Files, 2);
cycleStruct.cycleNS5FileInd = NaN(1,1); % which NS5 file this cycle belongs to
cycleStruct.cycleIndexWithinNS5File = NaN(1,1); % how the cycles are numbered in their original file
cycleStruct.cycleMinTime = NaN(1,1);
cycleStruct.cycleMaxTime = NaN(1,1);
cycleStruct.cycleMinAmp = NaN(1,1);
cycleStruct.cycleMaxAmp = NaN(1,1);
% MinMax & MaxMin Amp (these are used to compute slope, and can also be used to define a cycle's amplitude)
cycleStruct.cycleMinMaxAmp = NaN(1,1);
cycleStruct.cycleMaxMinAmp = NaN(1,1); 
cycleStruct.cycleMaxMaxStartMaxMinAmp = NaN(1,1); % REMEMBER: Start MaxMinAmp of MaxMaxCycle(N) = MaxMinAmp(N-1)
cycleStruct.cycleMinMinStartMinMaxAmp = NaN(1,1); % REMEMBER: Start MinMaxAmp of MinMinCycle(N) = MinMaxAmp(N)
% Duration
cycleStruct.cycleMinMinDuration = NaN(1,1);
cycleStruct.cycleMaxMaxDuration = NaN(1,1);
cycleStruct.cycleMinMaxDuration = NaN(1,1);
cycleStruct.cycleMaxMinDuration = NaN(1,1);
% Slopes: IMPORTANT IMPORTANT IMPORTANT
cycleStruct.cycleMinMaxSlope = NaN(1,1);
cycleStruct.cycleMaxMinSlope = NaN(1,1);
cycleStruct.cycleMaxMaxStartSlope = NaN(1,1); %   REMEMBER: StartSlope of MaxMaxCycle(N) = MaxMinSlope(N-1)
cycleStruct.cycleMinMinStartSlope = NaN(1,1); %   REMEMBER: StartSlope of MinMinCycle(N) = MinMaxSlope(N)
% Spectral Window that this cycle belongs to
cycleSpecStruct.cycleSpecWindowAssignWithinNS5File = NaN(1,1); % which spectral window this cycle belongs to
cycleSpecStruct.cycleSpecWindowAssign = NaN(1,1); % which spectral window this cycle belongs to
cycleSpecStruct.cycleSpecWindowBandPower = NaN(1,1);  % the power in the same band as this cycle

%% cycleSpikeStruct - stored with cell/lfp information
% Spike Count in each Cycle (from this cell)
cycleSpikeStruct.cycleMaxMinCount = NaN(1,1);
cycleSpikeStruct.cycleMinMaxCount = NaN(1,1);
cycleSpikeStruct.cycleMaxMaxCount = NaN(1,1);
cycleSpikeStruct.cycleMinMinCount = NaN(1,1);

%% spikePhaseStruct
spikePhaseStruct.spikeAssignMaxWithinNS5File  = NaN(1, numSpikes);
spikePhaseStruct.spikeAssignMinWithinNS5File  = NaN(1, numSpikes);
spikePhaseStruct.spikeAssignZeroWithinNS5File = NaN(1, numSpikes);
spikePhaseStruct.spikeAssignMax  = NaN(1, numSpikes);
spikePhaseStruct.spikeAssignMin  = NaN(1, numSpikes);
spikePhaseStruct.spikeAssignZero = NaN(1, numSpikes);
spikePhaseStruct.spikePhase = NaN(1, numSpikes);
spikePhaseStruct.spikeOffsets = NaN(1, numSpikes);
spikePhaseStruct.spikeOffsetsEnd = NaN(1, numSpikes);
spikePhaseStruct.spikeWithinCycleInd = NaN(1, numSpikes);

%% spikeRateCycleStruct - use for spike->cycle lookup. Length = length(spikeTimes)
spikeRateCycleStruct.spikeCycleMinTime  = NaN(1,numSpikes);
spikeRateCycleStruct.spikeCycleMaxTime  = NaN(1,numSpikes);
spikeRateCycleStruct.spikeCycleMinAmp    = NaN(1,numSpikes);
spikeRateCycleStruct.spikeCycleMaxAmp    = NaN(1,numSpikes);
% MinMax & MaxMin Amp (these are used to compute slope, and can also be used to define a cycle's amplitude)
spikeRateCycleStruct.spikeCycleMinMaxAmp = NaN(1,numSpikes); % remember: first half cycle is a minMax, as stored in humanProcessLFPMinMax
spikeRateCycleStruct.spikeCycleMaxMinAmp = NaN(1,numSpikes); % remember: second half cycle is a maxMin, as stored in humanProcessLFPMinMax
spikeRateCycleStruct.spikeCycleMaxMaxStartMaxMinAmp = NaN(1,numSpikes); % REMEMBER: Start MaxMinAmp of MaxMaxCycle(N) = MaxMinAmp(N-1)
spikeRateCycleStruct.spikeCycleMinMinStartMinMaxAmp = NaN(1,numSpikes); % REMEMBER: Start MinMaxAmp of MinMinCycle(N) = MinMaxAmp(N)
% Duration
spikeRateCycleStruct.spikeCycleMinMaxDuration = NaN(1,numSpikes);
spikeRateCycleStruct.spikeCycleMaxMinDuration = NaN(1,numSpikes);
spikeRateCycleStruct.spikeCycleMinMinDuration = NaN(1,numSpikes);
spikeRateCycleStruct.spikeCycleMaxMaxDuration = NaN(1,numSpikes);
% Slopes: IMPORTANT IMPORTANT IMPORTANT
spikeRateCycleStruct.spikeCycleMinMaxSlope = NaN(1,numSpikes);
spikeRateCycleStruct.spikeCycleMaxMinSlope = NaN(1,numSpikes);
spikeRateCycleStruct.spikeCycleMaxMaxStartSlope = NaN(1,numSpikes); %   REMEMBER: StartSlope of MaxMaxCycle(N) = MaxMinSlope(N-1)
spikeRateCycleStruct.spikeCycleMinMinStartSlope = NaN(1,numSpikes); %   REMEMBER: StartSlope of MinMinCycle(N) = MinMaxSlope(N)
% Spike Count in each Cycle (Rate is trivial to calculate so I will leave it out - just calculate it in the next level function when loading this data)
spikeRateCycleStruct.spikeCycleMinMaxCount = NaN(1,numSpikes);
spikeRateCycleStruct.spikeCycleMaxMinCount = NaN(1,numSpikes);
spikeRateCycleStruct.spikeCycleMaxMaxCount = NaN(1,numSpikes);
spikeRateCycleStruct.spikeCycleMinMinCount = NaN(1,numSpikes);

%% specStruct
specStruct.numWindowsPerNS5File = NaN(1,dataInfo(currSes).numNS5Files);
specStruct.ns5WindowMap = NaN(dataInfo(currSes).numNS5Files, 2);
specStruct.windowNS5FileInd = NaN(1,1); % which NS5 file this cycle belongs to
specStruct.windowIndexWithinNS5File = NaN(1,1); % how the cycles are numbered in their original file
specStruct.windowCenter = NaN(1,1);
specStruct.windowBandPower = NaN(1,1);
specStruct.windowCount = NaN(1,1);
specStruct.windowRate = NaN(1,1);
% specStruct.windowRateZ = NaN(1,1); % leaving this out because if I want Zscores I should do it later so that it incldes all the spikes across the concatenated windows

%% fullSpectrum
fullSpectrumStruct.f = NaN(1,1);
fullSpectrumStruct.windowCenter = NaN(1,1); % calling this windowCenter to make it compatible with specStruct
fullSpectrumStruct.S = NaN(1,1);

%% spikeRateWindowStruct - use for spike->spectral power lookups. Length = length(spikeTimes)
spikeRateWindowStruct.spikeSpecWindowAssignWithinNS5File = NaN(1,numSpikes); % which spectral window this spike belongs to
spikeRateWindowStruct.spikeSpecWindowAssign = NaN(1,numSpikes); % which spectral window this spike belongs to
spikeRateWindowStruct.spikeSpecWindowBandPower = NaN(1,numSpikes);  % the power in the same band as this spike
spikeRateWindowStruct.spikeSpecWindowCount = NaN(1, numSpikes); % the number of spikes in the spectral window that also contains this spike
%% THIS FILE RUNS AFTER human_spikeLFP_relationship_central_preprocessing.m
%  That file has already processed all the concatenated or single LFP for
%  single cylcle and spectrum calculations and stored the data in 3 types
%  of files:
%  1. overallZ_minmax or singleFileZ_minmax contain all the single cycle
%     data for a given channel and given frequency band (saved and loaded by humanProcessLFPMinMax)
%  2. filtered_LFP_nonZscored contains the filtered LFP if needed for plotting
%     or any other purpose (saved by humanProcessLFPMinMax)
%  3. spectrum data (saved and loaded by omarSpikeLFPRatePower (prefix is
%     not human because this file can work on any kind of input data)

%% GOALS OF THIS FILE
%  1. Load the single cycle LFP info and calc spike phase of each spike and
%     save in a single file for the entire concatenated time period
%  2. Call omarSpikeLFPRatePower, this time passing in the spikes. This
%     will assign the spikes to each spectral window. omarSpikeLFPRatePower function will
%     return, but not save the spike-spectrum relationship data. 
%  3. In the near future, add the cycle-offset code here as well
%  4. Don't worry about segmentation into awake-sleep-etc. here. The segementation will be
%     handled in the post-processing functions such as the one detailed
%     below.

%% GOALS of the NEXT FILE AFTER THIS ONE - which will process 1 or 2 cells
%  1. Read the phase/spikeAssignToSingleCycle/spikeAssigntoSpectralWindow
%     information from the file saved by
%     human_spike_LFP_relationship_central. This should also include the
%     rate of spikes per single cycle and per spectral window
%  2. Plot a very nice publication quality spikePhase histogram for one or
%     two cells.
%  3. Plot for 1 cell, the spike phase as a function of cycle amplitude
%     relationship for all spikes across all time windows or within a specified time period
%  4. Plot for 1 cell, the spike phase as a function of spectral window
%     power for all spikes across all time windows or within a specified time period
%  5. Now do overlaid plots of spike-kappa vs cycle amp or power for 2
%     cells - E and I
%  6. Now do overlaid plots of the spike-kappa vs cycle amp for multiple
%     brain states

%% DEFAULTS
% CONCATENATED_FILE_OUTPUT_DEFAULT_DIR = '../TEMP_ANALYSIS_DEPOT/concatenated_ns5s_for_clustering/'; % /outputParentDir/outputFilename .Where to store the data
% NS5_PROCESSED_DEFAULT_DIR = '../TEMP_ANALYSIS_DEPOT/processed_data'; %   % /subject/ns5filename/. All data related to this
% CELLS_PROCESSED_DEFAULT_DIR = '../TEMP_ANALYSIS_DEPOT/sorted_units/concatenated'; % /subject/
% FIGURES_DEFAULT_DIR = '../TEMP_ANALYSIS_DEPOT/figures/'; %   % /subject/ns5filename/. All data related to this
CONCATENATED_FILE_OUTPUT_DEFAULT_DIR = 'K:/concatenated_ns5s_for_clustering/'; % /outputParentDir/outputFilename .Where to store the data
NS5_PROCESSED_DEFAULT_DIR = 'K:/processed_data'; %   % /subject/ns5filename/. All data related to this
CELLS_PROCESSED_DEFAULT_DIR = 'H:/sorted_units/concatenated'; % /subject/   % leaving this on G: as the units are not that large but are accessed often, so best to keep it fast
FIGURES_DEFAULT_DIR = 'K:/figures/'; %   % /subject/ns5filename/. All data related to this

%% CELL/LFP/SESSION SETTINGS
session_list_human_lfp_spike_relationship_sfn_onwards;
% session_list_human_seizures_FOR_PAPER;
ses = [1]%[1 2 3 6 9 11 12 13 14]; %1,2=MG29;3=MG49;6=MG63;9=MG67;11,12=RIHE1,13=T1,14=T2
cellChToProcess = [1 43 51]; % leave empty to process all cells in session. Note that ALL cells on a given channel HAVE to be processed together because of the way I am automatically assigning the suffixes
cellMinQuality = 0; % a threshold for the minimum quality when allocating the allCellToProcessInSes variable
lfpChToProcess = []; % which LFP channels to do the phase locking of each cell with. 0 means only process the local LFP. If empty, it will read the channels from ns5VisualChs
maxElecDist = 400; % lfpChToProcess MUST be empty for this distance to be used. it is in microns. it will process all electrodes less than or equal to maxElecDist from the currCell's electrode
phaseMode = 2; % 1=min2minCycle, 2=max2maxCycle, 3=min2minCycle w/ asym, 4=max2maxCycle w/ asym
cellLFPPhaseInfoOverwrite = 0; % this will overwrite the main phaseInfo file if set to 1
cycleFileOverwrite = 0;
spikeSpecRateOverwrite = 0; % this is passed to omarSpikeLFPRatePower
saveFullSpectrum = 0;
fullSpectrumFileOverwrite = 0;
%% Settings for the spectrum
preFilterLFP = 1; % This will filter the LFP between 0.5 and min([900 Fs]) Hz before calculating the spectrogram
windowSize = 1; % NOTE: using 2 second window to allow 0.5 Hz to be better represented as well
tapers = [3 5]; %3 5 seems pointless, best to switch to 1 1?
fpass = [0.5 250];
numLogBins = 50;
%% Frequencies to process
%OLD ONE: lfpFreq = [0.5 1; 1 2; 2 3; 3 4; 4 5; 5 6; 6 7; 7 8; 8 9; 9 10; 10 11; 11 12; 12 14; 14 16; 16 18; 18 20; 20 22; 22 24; 24 26; 26 28; 28 30; 30 32; 32 34; 34 36; 36 38; 38 40; 40 42; 42 44; 44 46; 46 48; 48 50; 50 52; 52 54; 54 56; 56 58; 58 60; 60 62; 62 64; 64 66; 66 68; 68 70;];%[0.5 4; 4 8; 8 12; 12 18; 18 25; 25 32; 32 55; 65 110]; %[0.5 4; 4 8; 8 12; 12 18; 18 25; 25 30; 30 55; 65 110];
% NEW ONE:
lfpFreq = [150 250; 250 600]%[0.5 4; 4 8; 8 12; 30 55; 40 80; 0.5 1; 1 2; 2 3; 3 4; 4 5; 5 6; 6 7; 7 8; 8 9; 9 10; 10 11; 11 12; 12 14; 14 16; 16 18; 18 20; 20 22; 22 24; 24 26; 26 28; 28 30; 30 34; 34 38; 38 42; 42 46; 46 50; 50 54; 54 58; 58 62; 62 66; 66 70;]% 70 74; 74 78; 78 82; 82 86; 86 90];%[0.5 4; 4 8; 8 12; 12 18; 18 25; 25 32; 32 55; 65 110]; %[0.5 4; 4 8; 8 12; 12 18; 18 25; 25 30; 30 55; 65 110];
lfpFreqName = {'Delta', 'Theta', 'Alpha', 'Low Beta', 'High Beta', 'Beta-Gamma-Boundary', 'Low Gamma', 'High Gamma'};
plotLFPInd = [1 3 7];
plotLFPColorMap = [0.5 0.5 0.5; 0.8 0.2 0.2; 0.2 0.8 0.2];
numFreq = size(lfpFreq,1);

% Loop order:
% cell -> lfp -> filter bands -> single NS5s
% next question: store all spike info in one file or not?
% Yes, it makes the most sense. For each spike just keep spikeNS5FileInd
% and also keep a Nx2 matrix that keeps the firstSpikeInd and lastSpikeInd
% for each NS5 file. This will allow me to quickly go back and forth
% between the structures.
% For each spike, I want:
% - spikeTime, spikePhase, spikeAssignMin, spikeAssignMax, spikeCycleAmp (this is easy to compute from the minmax+spikeAssignMin),
%   spike

%% LOOP AROUND SESSIONS
for currSes = ses
    [arrayMap electrodeXY electrodeImp] = neuroportArrayData(dataInfo(currSes).subject);
    %% LOOP AROUND CELLS
    if isempty(cellChToProcess)
        cellIndToProcess = 1:length(dataInfo(currSes).cellChannel);
    else
        cellIndToProcess = find(ismember(dataInfo(currSes).cellChannel, cellChToProcess) == 1);
    end
    % also threshold for quality
    %cellIndToProcess = cellIndToProcess(find(dataInfo(currSes).cellQuality(cellIndToProcess) >= cellMinQuality)); % XXX THIS WILL SCREW UP THE a,b.
    prevCellCh = 0;
    numCellOnChan = 0;
    for currCellInd = cellIndToProcess
        currCellCh = dataInfo(currSes).cellChannel(currCellInd);
        if currCellCh ~= prevCellCh
            numCellOnChan = 1;
            spikeFile = [CELLS_PROCESSED_DEFAULT_DIR '/' dataInfo(currSes).subject '/' dataInfo(currSes).descriptiveFilename '/' 'chan' num2str(currCellCh) '/' dataInfo(currSes).descriptiveFilename '_ch' num2str(currCellCh) '.nex'];
            if ~exist(spikeFile, 'file')
                spikeFile = [CONCATENATED_FILE_OUTPUT_DEFAULT_DIR '/' dataInfo(currSes).subject '/' dataInfo(currSes).descriptiveFilename '/' dataInfo(currSes).descriptiveFilename '_ch' num2str(currCellCh) '.nex'];
                disp('USING spikeFile from the concatenated_ns5 directory');
            end
            info = ft_read_spike(spikeFile);
            FsSpike = info.hdr.FileHeader.Frequency;
            prevCellCh = currCellCh;
        else 
            numCellOnChan = numCellOnChan + 1;
        end
        currSuffix = convertNumberToLetter(numCellOnChan);
        chanLabel = ['Channel01' currSuffix];
        currCellName = sprintf('cell%g%s', currCellCh, currSuffix);
        withinChanCellInd = chan(info, chanLabel);
        spikeBins = double(info.timestamp{withinChanCellInd}); % XXX double is very important!!
        spikeTimes = spikeBins / FsSpike;
        numSpikes = length(spikeTimes);
        
        %% SPIKE TIMES ARE NOW READY - READ CONCAT METADATA so that spikes can be assigned to individual NS5 files
        % figure out the spike's NS5 metadata (for the concatenated NS5). do so by reading the metadata file        
        % XXX need to confirm that all concatenated LFP files have the
        % exact same number of bins, so that these spikes match up
        % perfectly with the LFPs from other channels. - I think this has to be the case
        concatNS5MetaDataFilename = sprintf('%s%s/%s/%s_ch%d_metatags.mat', CONCATENATED_FILE_OUTPUT_DEFAULT_DIR, dataInfo(currSes).subject, dataInfo(currSes).descriptiveFilename, dataInfo(currSes).descriptiveFilename, currCellCh);
        concatNS5MetaData = load(concatNS5MetaDataFilename);
        packetsPerFile = NaN(1, length(dataInfo(currSes).numNS5Files));
        FsPerFile = NaN(1, length(dataInfo(currSes).numNS5Files));
        for currNS5Ind = 1:dataInfo(currSes).numNS5Files % XXX Assuming this is the same as the contents of the metadata file
            packetsPerFile(currNS5Ind) = concatNS5MetaData.allMetaTags{currNS5Ind}.NumofPackets;
            FsPerFile(currNS5Ind) = concatNS5MetaData.allMetaTags{currNS5Ind}.SamplingFreq;
        end
        maxBinPerFile = cumsum(packetsPerFile); % now this can be used to find the spikeBins within the correct range
        maxTimePerFile = cumsum(packetsPerFile ./ FsPerFile);      
        
        %% LOOP AROUND LFP CHANNELS
        %  if LFP = 0, use only local spike channel, if empty use dataInfo.ns5VisChannels
        if isempty(lfpChToProcess), 
            if isempty(maxElecDist)
                allLFPChToProcess = dataInfo(currSes).ns5VisChannels;
            else
                allLFPChToProcess = sort(neuroportArrayGetElecAtDist(currCellCh, 0, maxElecDist, arrayMap, electrodeXY));
                % and intersect this set with the good ns5VisChannels to make sure bad channels are ignored
                allLFPChToProcess = intersect(allLFPChToProcess, dataInfo(currSes).ns5VisChannels);
                disp([currCellName ', MaxDist=' num2str(maxElecDist) ', LFPChToProcess=' num2str(allLFPChToProcess)]);
            end
        elseif lfpChToProcess == 0, allLFPChToProcess = currCellCh;
        else allLFPChToProcess = lfpChToProcess; end
        for currLFPCh = allLFPChToProcess
            %% LOOP AROUND FILTER BANDS
            for currFilterInd = 1:numFreq
                tic
                currLowFreq = lfpFreq(currFilterInd,1);
                currHighFreq = lfpFreq(currFilterInd,2);
                % spectrumLFPPostFix is defined in specData which is returend by omarSpikeLFPRatePower... I want to check for the existence of the phaseInfo
                % file BEFORE calling omarSpikeLFPRatePower so that I can skip over already existing ones. so it is important for me to redefine
                % spectrumLFPPostFix here...
                if preFilterLFP == 1, filterString = ['_prefiltered0.5-900'];
                else filterString = ''; end
                spectrumLFPPostfix = ['_window' num2str(windowSize*1000) 'ms_tapers' num2str(tapers(1)) '-' num2str(tapers(2)) '_fpass' num2str(fpass(1)) '-' num2str(fpass(2)) filterString];
                % NAMES OF THE SAVE/LOAD FILES
                saveDir = [NS5_PROCESSED_DEFAULT_DIR '/' dataInfo(currSes).subject '/' dataInfo(currSes).descriptiveFilename];
                cycleSaveFilename = sprintf('%s/lfp%g_freq%g-%g_mode%g_cycleinfo.mat', saveDir, currLFPCh, currLowFreq, currHighFreq, phaseMode);
                cycleSpecSaveFilename = sprintf('%s/lfp%g_freq%g-%g_mode%g%s_cyclespecinfo.mat', saveDir, currLFPCh, currLowFreq, currHighFreq, phaseMode, spectrumLFPPostfix);
                cellLFPSaveFilename = sprintf('%s/lfp%g_%s_freq%g-%g_mode%g%s_phaseinfo.mat', saveDir, currLFPCh, currCellName, currLowFreq, currHighFreq, phaseMode, spectrumLFPPostfix);
                
                %% EXISTING PHASEINFO FILE CHECK
                %  If the phaseInfo file exists, that means there is nothing to do here... just move on to the next LFP/Cell Pair
                if exist(cellLFPSaveFilename, 'file') & cellLFPPhaseInfoOverwrite == 0
                   continue;
                end
                
                %% LOOP AROUND EACH FILE
                % Make things cleaner by saving data in spikeNS5Struct, spikePhaseStruct, spikeRateCycleStruct, spikeRateSpecWindowStruct
                % To shorten this file, initiation done in helper .m file:
                human_spike_LFP_relationship_central_initialize_structs;
                spikeNS5Struct.spikeTimes = spikeTimes;
                spikeNS5Struct.numSpikes = numSpikes;
                numCyclesSoFar = 0;
                numWindowsSoFar = 0;
                
                %% EXISTING CYCLEINFO CHECK - loading the cycle file actually takes more time than just reassigning it, so edit this part out
                % if exist(cycleSaveFilename, 'file') & cycleFileOverwrite == 0
                %    load(cycleSaveFilename, 'cycleStruct'); % this will overwite the cycleStruct initialized above in human_..._initialize_structs
                %    reuseCycleStruct = 1;
                % else
                %    reuseCycleStruct = 0;
                % end
                
                %% LOOP AROUND ALL NS5 FILES
                for currNS5Ind = 1:dataInfo(currSes).numNS5Files % LOOP AROUND ALL NS5 FILES
                    %% PICK CURRENT SPIKE TIMES
                    currNS5Filename = dataInfo(currSes).ns5Filenames{currNS5Ind};
                    if currNS5Ind == 1
                        currSpikeInd = find(spikeBins > 0 & spikeBins < maxBinPerFile(currNS5Ind));
                        spikeZeroBin = 0; spikeZeroTime = 0;
                    else
                        currSpikeInd = find(spikeBins > maxBinPerFile(currNS5Ind-1) & spikeBins <= maxBinPerFile(currNS5Ind));
                        spikeZeroBin = maxBinPerFile(currNS5Ind-1);
                        spikeZeroTime = spikeZeroBin / FsSpike;
                    end
                    currSpikeTimes = spikeTimes(currSpikeInd) - spikeZeroTime;
                    currSpikeBins  = spikeBins(currSpikeInd) - spikeZeroBin;
                    % doublecheck to make we are never skipping spike inds (all the spikes in one NS5 are contiguous)
                    if max(unique(diff(currSpikeInd))) > 1 error(['NS5Ind ' currNS5Ind ' is returning non-contiguous spike indices somehow.']); end
                    %% STORE spikeNS5Struct related info
                    spikeNS5Struct.spikeNS5FileInd(currSpikeInd) = currNS5Ind;
                    if ~isempty(currSpikeInd)
                        spikeNS5Struct.ns5SpikeMap(currNS5Ind,1) = currSpikeInd(1);
                        spikeNS5Struct.ns5SpikeMap(currNS5Ind,2) = currSpikeInd(end);
                    end
                    spikeNS5Struct.spikeTimesWithinNS5File(currSpikeInd) = currSpikeTimes;
                    spikeNS5Struct.packetsPerFile = packetsPerFile;
                    spikeNS5Struct.FsPerFile = FsPerFile;
                    spikeNS5Struct.maxBinPerFile = maxBinPerFile;
                    spikeNS5Struct.maxTimePerFile = maxTimePerFile;
                    
                    %% SPIKE PHASE
                    % Correct spikes for this NS5Ind have been picked - now process and store them in "curr" variables before
                    % transferring to the permananet concatenated variables
                    outputString = sprintf('PHASE/RATE: About to Process Freq=%g-%g, File %s, NS5Ind=%g, Cell=%g%s, LFP=%g', currLowFreq, currHighFreq, dataInfo(currSes).ns5Filenames{currNS5Ind}, currNS5Ind, currCellCh, currSuffix, currLFPCh); disp(outputString);
                    currProcessedDataDir = sprintf('%s/%s/%s/', NS5_PROCESSED_DEFAULT_DIR, dataInfo(currSes).subject, dataInfo(currSes).ns5Filenames{currNS5Ind});
                    lfpFilenamePrefix = [currProcessedDataDir 'lfp' num2str(currLFPCh) '_dec'];
                    useOverallMeanStd = 1;
                    cycleInfo = humanProcessLFPMinMax([], [], [], currLowFreq, currHighFreq, lfpFilenamePrefix, 1,0,0, useOverallMeanStd);
                    [currSpikeAssignMax currSpikeAssignMin currSpikeAssignZero] = assignSpikesToCycles(currSpikeTimes, ...
                        cycleInfo.ripMinTimes, cycleInfo.ripMaxTimes);
                    [currSpikePhase currSpikeOffsets currSpikeOffsetsEnd] = assignPhaseToSpikes(currSpikeTimes, currSpikeAssignMax, currSpikeAssignMin, ...
                        cycleInfo.ripMinTimes, cycleInfo.ripMaxTimes, phaseMode, 1, currSpikeAssignZero);
                    [currWithinCycleInd checkNumSpikesInCycle] = assignWithinCycleIndexToSpikes(currSpikePhase, currSpikeAssignMax, currSpikeAssignMin, phaseMode);
                    % checkNumSpikesInCycle is the same as spikeRateCycleStruct.spikeCycleMaxMaxCount and can ignored or used for debugging as necessary.
                    %% SINGLE CYCLE RATE
                    currSpikeRateStruct = assignRatesToCycles(currSpikeTimes, currSpikeAssignMax, currSpikeAssignMin, cycleInfo.ripMinTimes, cycleInfo.ripMaxTimes);
                    
                    %% SPECTRUM AND SPECTRAL WINDOW RATE
                    [specData spikeSpecData] = omarSpikeLFPRatePower(currSpikeTimes, FsSpike, currCellName, [], [], [], lfpFilenamePrefix, windowSize, tapers, fpass, numLogBins, preFilterLFP, spikeSpecRateOverwrite);
                    % make sure the spectrumLFPPostfix returned by omarSpikeLFPRatePower is the same as that defined above in this function
                    if strcmp(specData.spectrumLFPPostfix, spectrumLFPPostfix) == 0, error('MISMATCHING spectrumLFPPostfix'); end;
                    fInd = find(specData.f >= currLowFreq & specData.f <= currHighFreq);
                    windowBandPower = nanmean(specData.S(:,fInd),2);
                                     
                    %% STORE CYCLE INFO in cycleStruct
                    cycleStruct.numCyclesPerNS5File(currNS5Ind) = length(cycleInfo.ripMinTimes);
                    cycleStruct.ns5CycleMap(currNS5Ind,1) = numCyclesSoFar + 1;
                    cycleStruct.ns5CycleMap(currNS5Ind,2) = numCyclesSoFar + length(cycleInfo.ripMinTimes);
                    currCycleInd = cycleStruct.ns5CycleMap(currNS5Ind,1):cycleStruct.ns5CycleMap(currNS5Ind,2);
                    cycleStruct.cycleNS5FileInd(currCycleInd) = currNS5Ind;
                    cycleStruct.cycleIndexWithinNS5File(currCycleInd) = 1:length(cycleInfo.ripMinTimes);
                    cycleStruct.cycleMinTime(currCycleInd) = cycleInfo.ripMinTimes + spikeZeroTime;
                    cycleStruct.cycleMaxTime(currCycleInd) = cycleInfo.ripMaxTimes + spikeZeroTime;
                    cycleStruct.cycleMinAmp(currCycleInd)  = cycleInfo.ripMinAmp;
                    cycleStruct.cycleMaxAmp(currCycleInd)  = cycleInfo.ripMaxAmp;
                    % MinMax & MaxMin Amp (these are used to compute slope, and can also be used to define a cycle's amplitude)
                    cycleStruct.cycleMinMaxAmp(currCycleInd) = cycleInfo.ripMaxAmp - cycleInfo.ripMinAmp;
                    cycleStruct.cycleMaxMinAmp(currCycleInd(1:end-1)) = cycleInfo.ripMaxAmp(1:end-1) - cycleInfo.ripMinAmp(2:end); % last one is NaN
                    cycleStruct.cycleMaxMinAmp(currCycleInd(end)) = NaN; % explicitly setting the last bin to NaN to make it clear that this has 1 less bin
                    cycleStruct.cycleMinMinStartMinMaxAmp(currCycleInd) = cycleStruct.cycleMinMaxAmp(currCycleInd); % REMEMBER: Start MinMaxAmp of MinMinCycle(N) = MinMaxAmp(N)
                    cycleStruct.cycleMaxMaxStartMaxMinAmp(currCycleInd(2:end)) = cycleStruct.cycleMinMinStartMinMaxAmp(currCycleInd(1:end-1));% REMEMBER: Start MaxMinAmp of MaxMaxCycle(N) = MaxMinAmp(N-1)
                    cycleStruct.cycleMaxMaxStartMaxMinAmp(currCycleInd(1)) = NaN; % explicitly set to NaN - otherwise it will be stored as 0 when Matlab extends the array, since this array grows in every loop
                    % Duration
                    cycleStruct.cycleMinMinDuration(currCycleInd) = cycleInfo.ripDurations;
                    cycleStruct.cycleMaxMaxDuration(currCycleInd) = cycleInfo.ripMaxMaxDurations;
                    cycleStruct.cycleMinMaxDuration(currCycleInd) = cycleInfo.ripMinMaxDurations;
                    cycleStruct.cycleMaxMinDuration(currCycleInd(1:end-1)) = cycleInfo.ripMaxMinDurations(1:end); % ripMaxMinDurations is 1 bin shorter than all the other cycle variables because there is no end MaxMinCycle
                    cycleStruct.cycleMaxMinDuration(currCycleInd(end)) = NaN; % explicitly setting the last bin to NaN to make it clear that this has 1 less bin
                    % Slopes: IMPORTANT IMPORTANT IMPORTANT
                    cycleStruct.cycleMinMaxSlope(currCycleInd) = cycleStruct.cycleMinMaxAmp(currCycleInd) ./ cycleStruct.cycleMinMaxDuration(currCycleInd);
                    cycleStruct.cycleMaxMinSlope(currCycleInd) = cycleStruct.cycleMaxMinAmp(currCycleInd) ./ cycleStruct.cycleMaxMinDuration(currCycleInd);
                    cycleStruct.cycleMinMinStartSlope(currCycleInd) = cycleStruct.cycleMinMaxSlope(currCycleInd); % REMEMBER: StartSlope of MinMinCycle(N) = MinMaxSlope(N)
                    cycleStruct.cycleMaxMaxStartSlope(currCycleInd(2:end)) = cycleStruct.cycleMaxMinSlope(currCycleInd(1:end-1)); % REMEMBER: StartSlope of MaxMaxCycle(N) = MaxMinSlope(N-1)
                    cycleStruct.cycleMaxMaxStartSlope(currCycleInd(1)) = NaN; % explicitly set to NaN - otherwise it will be stored as 0 when Matlab extends the array, since this array grows in every loop
                    
                    % Spectral Window that this cycle belongs to
                    % Assign the minima and maxima of each cycle to a spectral time window
                    [cycleWindowAssign ignoreSpikeAssignMin] = assignSpikesToCycles(cycleInfo.ripMinTimes, specData.time.timeBinEdges(1:end-1), specData.time.timeBinCenters); % treat the time windows as min-to-min cycles with timeBinCenters being the max time in the middle.... so for a min-2-min cycle we only care about spikeAssignMax
                    currCycleNonNan = find(isnan(cycleWindowAssign) == 0);
                    currCycleNan = find(isnan(cycleWindowAssign) == 1);
                    cycleSpecStruct.cycleSpecWindowAssignWithinNS5File(currCycleInd) = cycleWindowAssign;
                    cycleSpecStruct.cycleSpecWindowAssign(currCycleInd) = cycleWindowAssign + numWindowsSoFar;
                    cycleSpecStruct.cycleSpecWindowBandPower(currCycleInd(currCycleNonNan)) = windowBandPower(cycleWindowAssign(currCycleNonNan));
                    cycleSpecStruct.cycleSpecWindowBandPower(currCycleInd(currCycleNan)) = NaN;

                    
                    %% Spike Count in each Cycle (from this cell)
                    cycleSpikeStruct.cycleMaxMinCount(currCycleInd) = currSpikeRateStruct.maxMinCount;
                    cycleSpikeStruct.cycleMinMaxCount(currCycleInd) = currSpikeRateStruct.minMaxCount;
                    cycleSpikeStruct.cycleMaxMaxCount(currCycleInd) = currSpikeRateStruct.maxMaxCount;
                    cycleSpikeStruct.cycleMinMinCount(currCycleInd) = currSpikeRateStruct.minMinCount;
                    
                    %% STORE SPIKE PHASE in the OVERALL CONCAT vectors/struct
                    spikePhaseStruct.spikeAssignMaxWithinNS5File(currSpikeInd)  = currSpikeAssignMax;
                    spikePhaseStruct.spikeAssignMinWithinNS5File(currSpikeInd)  = currSpikeAssignMin;
                    spikePhaseStruct.spikeAssignZeroWithinNS5File(currSpikeInd) = currSpikeAssignZero;
                    spikePhaseStruct.spikeAssignMax(currSpikeInd)  = currSpikeAssignMax + numCyclesSoFar;
                    spikePhaseStruct.spikeAssignMin(currSpikeInd)  = currSpikeAssignMin + numCyclesSoFar;
                    spikePhaseStruct.spikeAssignZero(currSpikeInd) = currSpikeAssignZero + numCyclesSoFar;
                    spikePhaseStruct.spikePhase(currSpikeInd)      = currSpikePhase;
                    spikePhaseStruct.spikeOffsets(currSpikeInd)    = currSpikeOffsets;
                    spikePhaseStruct.spikeOffsetsEnd(currSpikeInd) = currSpikeOffsetsEnd; % something is wrong with the way offsetsEnd is being calculated in assignPhaseToSpikes - XXX look into it
                    spikePhaseStruct.spikeWithinCycleInd(currSpikeInd) = currWithinCycleInd;
                                        
                    %% STORE SPIKE->CYCLE mapping info in spikeRateCycleStruct
                    currSpikeNonNan = find(isnan(currSpikePhase) == 0 & isnan(currSpikeAssignMax) == 0); % The currSpikeAssignMax check is for a very special case where a spike is in the last cycle and phaseMode=2 (i.e. phase is calculated max-2-max). This spike will have a phase, but it might belong to the very last cycleMax - which I do not allow. Hence it's spikeAssignMax will be NaN (see assignSpikesToCycle). Including this check is a nice little precaution and prevents an error in the lines below.
                    % Save the amplitude/slope/duration/rate of the cycle in which each spike occurred - this makes
                    % reverse lookup (spike->cycle property) much quicker and easier 
                    spikeRateCycleStruct.spikeCycleMinTime(currSpikeInd(currSpikeNonNan)) = cycleInfo.ripMinTimes(currSpikeAssignMin(currSpikeNonNan));
                    spikeRateCycleStruct.spikeCycleMaxTime(currSpikeInd(currSpikeNonNan)) = cycleInfo.ripMaxTimes(currSpikeAssignMax(currSpikeNonNan));
                    spikeRateCycleStruct.spikeCycleMinAmp(currSpikeInd(currSpikeNonNan)) = cycleInfo.ripMinAmp(currSpikeAssignMin(currSpikeNonNan));
                    spikeRateCycleStruct.spikeCycleMaxAmp(currSpikeInd(currSpikeNonNan)) = cycleInfo.ripMaxAmp(currSpikeAssignMax(currSpikeNonNan));
                    % MinMax & MaxMin Amp (these are used to compute slope, and can also be used to define a cycle's amplitude)
                    spikeRateCycleStruct.spikeCycleMinMaxAmp(currSpikeInd(currSpikeNonNan)) = cycleStruct.cycleMinMaxAmp(spikePhaseStruct.spikeAssignMax(currSpikeInd(currSpikeNonNan)));
                    spikeRateCycleStruct.spikeCycleMaxMinAmp(currSpikeInd(currSpikeNonNan)) = cycleStruct.cycleMaxMinAmp(spikePhaseStruct.spikeAssignMax(currSpikeInd(currSpikeNonNan)));
                    spikeRateCycleStruct.spikeCycleMinMinStartMinMaxAmp(currSpikeInd(currSpikeNonNan)) = cycleStruct.cycleMinMinStartMinMaxAmp(spikePhaseStruct.spikeAssignMax(currSpikeInd(currSpikeNonNan))); % using spikeAssignMax since this is a minMin cycle
                    spikeRateCycleStruct.spikeCycleMaxMaxStartMaxMinAmp(currSpikeInd(currSpikeNonNan)) = cycleStruct.cycleMaxMaxStartMaxMinAmp(spikePhaseStruct.spikeAssignMin(currSpikeInd(currSpikeNonNan))); % using spikeAssignMin since this is a maxMax cycle
                    % Duration
                    spikeRateCycleStruct.spikeCycleMinMaxDuration(currSpikeInd(currSpikeNonNan)) = cycleInfo.ripMinMaxDurations(currSpikeAssignMax(currSpikeNonNan));
                    spikeRateCycleStruct.spikeCycleMaxMinDuration(currSpikeInd(currSpikeNonNan)) = cycleInfo.ripMaxMinDurations(currSpikeAssignMax(currSpikeNonNan));
                    spikeRateCycleStruct.spikeCycleMinMinDuration(currSpikeInd(currSpikeNonNan)) = cycleInfo.ripDurations(currSpikeAssignMax(currSpikeNonNan));
                    spikeRateCycleStruct.spikeCycleMaxMaxDuration(currSpikeInd(currSpikeNonNan)) = cycleInfo.ripMaxMaxDurations(currSpikeAssignMin(currSpikeNonNan));
                    % Slope
                    spikeRateCycleStruct.spikeCycleMinMaxSlope(currSpikeInd(currSpikeNonNan)) = cycleStruct.cycleMinMaxSlope(spikePhaseStruct.spikeAssignMax(currSpikeInd(currSpikeNonNan)));
                    spikeRateCycleStruct.spikeCycleMaxMinSlope(currSpikeInd(currSpikeNonNan)) = cycleStruct.cycleMaxMinSlope(spikePhaseStruct.spikeAssignMax(currSpikeInd(currSpikeNonNan)));
                    spikeRateCycleStruct.spikeCycleMinMinStartSlope(currSpikeInd(currSpikeNonNan)) = cycleStruct.cycleMinMinStartSlope(spikePhaseStruct.spikeAssignMax(currSpikeInd(currSpikeNonNan)));
                    spikeRateCycleStruct.spikeCycleMaxMaxStartSlope(currSpikeInd(currSpikeNonNan)) = cycleStruct.cycleMaxMaxStartSlope(spikePhaseStruct.spikeAssignMin(currSpikeInd(currSpikeNonNan))); % using spikeAssignMin since this is for a maxMax cycle
                    % Count
                    spikeRateCycleStruct.spikeCycleMinMaxCount(currSpikeInd(currSpikeNonNan)) = cycleSpikeStruct.cycleMinMaxCount(spikePhaseStruct.spikeAssignMax(currSpikeInd(currSpikeNonNan)));
                    spikeRateCycleStruct.spikeCycleMaxMinCount(currSpikeInd(currSpikeNonNan)) = cycleSpikeStruct.cycleMaxMinCount(spikePhaseStruct.spikeAssignMax(currSpikeInd(currSpikeNonNan)));
                    spikeRateCycleStruct.spikeCycleMinMinCount(currSpikeInd(currSpikeNonNan)) = cycleSpikeStruct.cycleMinMinCount(spikePhaseStruct.spikeAssignMax(currSpikeInd(currSpikeNonNan)));
                    spikeRateCycleStruct.spikeCycleMaxMaxCount(currSpikeInd(currSpikeNonNan)) = cycleSpikeStruct.cycleMaxMaxCount(spikePhaseStruct.spikeAssignMin(currSpikeInd(currSpikeNonNan))); % using spikeAssignMin since this is for a maxMax cycle
                    
                    %% STORE SPECTRAL POWER AND RELATED STUFF
                    % calculate the mean power in the band of interest
                    specStruct.numWindowsPerNS5File(currNS5Ind) = length(specData.t);
                    specStruct.ns5WindowMap(currNS5Ind,1) = numWindowsSoFar + 1;
                    specStruct.ns5WindowMap(currNS5Ind,2) = numWindowsSoFar + length(specData.t);
                    currWindowInd = specStruct.ns5WindowMap(currNS5Ind,1):specStruct.ns5WindowMap(currNS5Ind,2);
                    specStruct.windowNS5FileInd(currWindowInd) = currNS5Ind;
                    specStruct.windowIndexWithinNS5File(currWindowInd) = 1:length(specData.t);
                    specStruct.windowCenter(currWindowInd) = specData.t + spikeZeroTime;
                    specStruct.windowBandPower(currWindowInd) = windowBandPower;
                    specStruct.windowCount(currWindowInd) = spikeSpecData.timeBinSpikeCount;
                    specStruct.windowRate(currWindowInd) = spikeSpecData.timeBinRate;                    
                    
                    %% STORE FULL SPECTRUM
                    if saveFullSpectrum == 1
                        fullSpectrumStruct.f = specData.f; % this will just get overriden for each ns5Ind
                        fullSpectrumStruct.windowCenter(currWindowInd) = specData.t + spikeZeroTime;
                        fullSpectrumStruct.S(1:length(specData.f), currWindowInd) = specData.S'; % stored spectrum is fxt (freq in rows, time in columns)
                    end
                    
                    %% STORE SPIKE->SPEC WINDOW mapping info in spikeRateWindowStruct
                    currSpikeWindowAssignNonNan = find(isnan(spikeSpecData.spikeTimeAssignMax) == 0);
                    spikeRateWindowStruct.spikeSpecWindowAssignWithinNS5File(currSpikeInd) = spikeSpecData.spikeTimeAssignMax;
                    spikeRateWindowStruct.spikeSpecWindowAssign(currSpikeInd) = spikeSpecData.spikeTimeAssignMax + numWindowsSoFar;
                    spikeRateWindowStruct.spikeSpecWindowBandPower(currSpikeInd(currSpikeWindowAssignNonNan)) = windowBandPower(spikeSpecData.spikeTimeAssignMax(currSpikeWindowAssignNonNan));
                    spikeRateWindowStruct.spikeSpecWindowCount(currSpikeInd(currSpikeWindowAssignNonNan)) = specStruct.windowCount(spikeRateWindowStruct.spikeSpecWindowAssign(currSpikeInd(currSpikeWindowAssignNonNan)));
                                        
                    %% INCREMENT THE CYCLES AND WINDOWS
                    numCyclesSoFar = numCyclesSoFar + length(cycleInfo.ripMinTimes);
                    numWindowsSoFar = numWindowsSoFar + length(specData.t);
                end % End loop around NS5 files that have been concatenated together
                                
                %% Save this file in the processed_data/subject/descriptiveFilename/ directory
                if ~exist(saveDir, 'dir')
                    mkdir(saveDir);
                end
                if ~exist(cycleSaveFilename, 'file') | (cycleFileOverwrite == 1 & numCellOnChan == 1) % if overwriting, only overwrite for the first cell on each channel
                    save(cycleSaveFilename, 'cycleStruct', 'currNS5Ind');
                end
                if ~exist(cycleSpecSaveFilename, 'file') | (cycleFileOverwrite == 1 & numCellOnChan == 1) % if overwriting, only overwrite for the first cell on each channel
                    save(cycleSpecSaveFilename, 'cycleSpecStruct', 'currNS5Ind');
                end
                % SAVE fhe PHASEINFO - the overwrite checks have already been performed above.
                save(cellLFPSaveFilename, 'spikeNS5Struct', 'spikePhaseStruct', 'spikeRateCycleStruct', 'spikeRateWindowStruct', 'cycleSpikeStruct', 'specStruct', 'currNS5Ind');
                
                % Now save the full spectrum if that option has been selected
                if saveFullSpectrum == 1
                    fullSpectrumSaveFilename = sprintf('%s/lfp%g%s_fullspectruminfo.mat', saveDir, currLFPCh, specData.spectrumLFPPostfix);
                    if ~exist(fullSpectrumSaveFilename, 'file') | (fullSpectrumFileOverwrite == 1 & numCellOnChan == 1) % if overwriting, only overwrite for the first cell on each channel
                        save(fullSpectrumSaveFilename, 'fullSpectrumStruct');
                    end
                end
                clear cycleStruct cycleSpecStruct spikeNS5Struct spikePhaseStruct spikeRateCycleStuct spikeRateWindowStruct cycleSpikeStruct specStruct
                toc
            end % End loop around filter bands
        end % End loop around LFP channels
    end % End loop around Cells
end % End loop around ses
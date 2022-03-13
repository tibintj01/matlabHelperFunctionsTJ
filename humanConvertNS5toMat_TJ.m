function humanConvertNS5toMat(checkBeforeConverting, subject, ns5Filename, ns5VisChannels, ns5SeqChannels, originalDataDir, processedDataDir, saveDecimated)
%
% function humanConvertNS5toMat(checkBeforeConverting, subject, ns5Filename, ns5VisChannels, ns5SeqChannels, originalDataDir, processedDataDir, decimate)
% 
% Omar Ahmed, Original date: August 15, 2011. Turned into a function on
% December 27, 2011. It can be called from human_GRAND_preprocessing_script
% or called from the humanLoadNS5Channels.m
% 
% Input Parameters (only subject, ns5Filename, ns5VisChannels are needed)
% -----------------------------------------------------------------------
% checkBeforeConverting: If 1, the code will check for the existence of the
%                        original and decimiated files and only convert if 
%                        the processed .mat files don't already exist.
%                        If 0, the conversion will go ahead no matter what.
% subject:  'MG49' or 'MG29' or 'BW11', etc.
% ns5Filename:  '20110615-094530-023', no need for the .ns5 file extension.
%               This function only processes .ns5 files. Ignore the following
%               as NEVs have flawed timestamps and are useless: Note that both the
%               NS5 file and the corresponding NEV file BOTH need to be
%               located in the same folder [The NEV is absolutely necessary
%               to extract timestamp information so that multiple NS5 files
%               can be strung together with the correct blank space in
%               between. This allows them to be correctly aligned with spikes
%               that were from multiple NEV files, but were clustered together
%               as one].
% ns5VisChannels: 0 means process all channels in file.
%                 Otherwise these channel numbers are the ones that appear on
%                 Blackrock Central. These are the channel numbers
%                 that will used to save the channel data so that
%                 from now on they can be loaded more easily. If the
%                 Blackrock Central channel numbers don't start at 1 and
%                 are not consecutive, then you also need to pass in the
%                 corresponding ns5SeqChannels. 
% ns5SeqChannels: [OPTIONAL] Defaults to the same values as the ns5VisChannels.
%                 Otherwise these are sequential channel numbers.
%                 E.g. if the channel numbers range from 33-64 on Blackrock
%                 Central, then the channels passed in should actually be
%                 1:32. These are the numbers that will be used to load the
%                 data from the NS5 {DOUBLE-CHECK that this is the correct
%                 thing to do). Checked - it is correct.
% originalDataDir: [OPTIONAL] If the NS5/NEV original files should be
%                  loaded from a custom location
% processedDataDir: [OPTIONAL] If the processed files should be stored in a
%                   custom location                   
% decimate: [OPTIONAL] Defaults to 1. Whether or not to also generate decimated version of the data.
%           The default factor for decimation is 15, resulting in Fs=2000

NS5_ORIGINAL_DEFAULT_DIR        = 'H:/AHMEDLAB_DATA_COPIES'; % /subject/NS5_NEV/. Be sure to keep this folder structure. Original NS5s and NEVs are both kept in this folder for this subject
NS5_ORIGINAL_DEFAULT_DIR_BACKUP = 'I:/AHMEDLAB_DATA_COPIES'; % Some patients are stored on I: instead of H:
NS5_PROCESSED_DEFAULT_DIR       = 'J:/processed_data';   % /subject/ns5filename/. All data related to this 
numSimChLoad = 100; % How many channels to process at a time, 5 is optimal to not run out of memory
decFactor = 15;   % Decimate by a factor of 15, resulting in a 2000 Hz Fs (original = 30000 for this dataset

%% SET DEFAULT VALUES FOR OPTIONAL PARAMETERS
if exist('ns5SeqChannels') ~= 1 || isempty(ns5SeqChannels)
    disp('No sequential NS5 channel numbers passed into this function. Defaulting to the visual NS5 channel numbers');
    ns5SeqChannels = ns5VisChannels; % if no sequential channel numbers are passed in, default to the same as the visual numbers, and assume that there is a 1:1 correspondence between the two
end

if exist('originalDataDir') ~= 1
    disp(['Using Default ORIGINAL Data Dir: ' NS5_ORIGINAL_DEFAULT_DIR]);
    disp(['Using Default ORIGINAL BACKUP Data Dir: ' NS5_ORIGINAL_DEFAULT_DIR_BACKUP]);
    originalDataDir = NS5_ORIGINAL_DEFAULT_DIR;
    originalDataDirBackup = NS5_ORIGINAL_DEFAULT_DIR_BACKUP;
end

if exist('processedDataDir') ~= 1
    disp(['Using Default PROCESSED Data Dir: ' NS5_PROCESSED_DEFAULT_DIR]);
    processedDataDir = NS5_PROCESSED_DEFAULT_DIR;
end

if exist('saveDecimated') ~= 1
    disp('Will decimate data');
    saveDecimated = 1;
end

%% DIRECTORY AND FILENAME SETUP
% dataDir = sprintf('../data/%s/', subject);
originalDataDir = sprintf('%s/%s/NS5_NEV/', originalDataDir, subject)
%originalDataDir = sprintf('%s/%s/NS5/', originalDataDir, subject)
%fds
fullPath = [originalDataDir ns5Filename '.ns5'];
if exist('originalDataDirBackup') == 1 % Only if no originalDataDir was passed into this function
    originalDataDirBackup = sprintf('%s/%s/NS5_NEV/', originalDataDirBackup, subject);
    fullPathBackup = [originalDataDirBackup ns5Filename '.ns5'];
end
processedDataDir = sprintf('%s/%s/%s/', processedDataDir, subject, ns5Filename);
if exist(processedDataDir) ~= 7
    mkdir(processedDataDir);
end

%% CHECK IF THE ORIGINALD_DATA_DIR or ORIGINAL_DATA_DIR_BACKUP contains the NS5 file before proceeding
if exist('originalDataDirBackup') == 1 % figure out if the NS5 file exists in the original or originalBackup directory
    if exist(fullPath) ~= 2 % if the file is not in the original data folder
        if exist(fullPathBackup) ~= 2
            error(['File not found in either ORIGINAL or ORIGINAL_BACKUP Data Dir. Exiting']);
        else % the file exists in the backup dir, so use that, and all is good
            disp(['Switching to Default ORIGINAL_BACKUP Data Dir: ' NS5_ORIGINAL_DEFAULT_DIR_BACKUP]);
            fullPath = fullPathBackup;
        end
    end
end

%% If ns5VisChannels == 0, then read the header to find out how many channels there are in total
if ns5VisChannels == 0
    allLFPData = openNSx('report', fullPath); % just read the header to get the full channel list
    ns5VisChannels =  allLFPData.MetaTags.ChannelID'; % this sets it to all the channels
    % This will return [1:96 129 130]. Therefore it is returning visual
    % channel numbers.
    ns5SeqChannels = 1:length(ns5VisChannels);
end

%% CHECK TO SEE IF processed .mat files already exist
%  channels is the variable that will be used to do the conversion from NS5
%  to original.mat files
channelIndToProcess = []; % these are indices into ns5SeqChannels and ns5VisChannels
if checkBeforeConverting == 1
    for i=1:length(ns5SeqChannels)
        currVisCh = ns5VisChannels(i); % all filenames are based on visual channel numbers
        saveFilename = ['lfp' num2str(currVisCh) '_original.mat'];
        saveFullPath = [processedDataDir saveFilename];
        if exist(saveFullPath, 'file') == 0
            % File Doesn't Exist. Add to "channels" list
            channelIndToProcess = [channelIndToProcess i]; 
        end
    end
end

%% Now process N channels at a time - 5 seems optimal so as to not run out of memory
%  Note that if you try to load channels [6 97], then openNSx will actually read ALL channels in between as well!!
%  This is pretty much a guaranteed way to crash matlab and your computer. So the calling functions should make sure
%  that horribly discontinuous channels are not passed into 
numLoadRounds = ceil(length(channelIndToProcess)/numSimChLoad);
for currRound=1:numLoadRounds
    if currRound == numLoadRounds
        currSeqChannels = ns5SeqChannels(channelIndToProcess((numSimChLoad*(currRound-1)+1):length(channelIndToProcess)))
        currVisChannels = ns5VisChannels(channelIndToProcess((numSimChLoad*(currRound-1)+1):length(channelIndToProcess)))
    else
        currSeqChannels = ns5SeqChannels(channelIndToProcess((numSimChLoad*(currRound-1)+1):numSimChLoad*currRound))
        currVisChannels = ns5VisChannels(channelIndToProcess((numSimChLoad*(currRound-1)+1):numSimChLoad*currRound))
    end
    disp(['Round ' num2str(currRound) '/' num2str(numLoadRounds) ': Loading Visual Channel(s) ' num2str(currVisChannels)]);
    allLFPData = openNSx('report', 'read', fullPath, ['c:' num2str(currSeqChannels)]); % Note the use of Seq Channels when calling openNSx
    metaTags = allLFPData.MetaTags;
    Fs = metaTags.SamplingFreq;

    % and save these N channels before moving onto reading the next round of N
    for i=1:length(currSeqChannels)
        currVisCh = currVisChannels(i);
        currSeqCh = currSeqChannels(i);
        lfp = allLFPData.Data(i,:);
        % Calcualate the LFP stats necessary to compute the composite standard deviation
        % See http://www.burtonsys.com/climate/composite_standard_deviations.html
        clear lfpStats;
        lfpStats.mean = nanmean(lfp);
        lfpStats.var = nanvar(double(lfp));
        lfpStats.std = sqrt(lfpStats.var);
        lfpStats.N = length(lfp);
        lfpStats.nonNanN = length(find(isnan(lfp) == 0));
        lfpStats.errSumSq = lfpStats.var * (lfpStats.N - 1); % this is important for calculating the overall standard deviation
        lfpStats.rms = sqrt(nanmean(lfp.^2));
        lfpStats.max = nanmax(lfp);
        lfpStats.min = nanmin(lfp);
        
        saveFilename = ['lfp' num2str(currVisCh) '_original.mat'];
        saveFullPath = [processedDataDir saveFilename];
        disp(['Saving LFP Visual Channel ' num2str(currVisCh)]); tic;
        visChNum = currVisCh;
        seqChNum = currSeqCh;
        numPackets = length(lfp);
        save(saveFullPath, 'lfp', 'Fs', 'metaTags', 'numPackets', 'visChNum', 'seqChNum', 'lfpStats'); % the decimated file should contain the same info, PLUS timeAxis (timeAxis not needed here as we are not skipping anything)
        toc;

        clear lfp
    end
    clear allLFPData Fs metaTags
end

%Tibin
return
if saveDecimated ~= 1
    return;
end

%% CHECK TO SEE IF processed and DECIMATED .mat files already exist
%  channels is the variable that will be used to do the conversion from NS5
%  to original.mat files
channelIndToProcess = []; % these are indices into ns5SeqChannels and ns5VisChannels
if checkBeforeConverting == 1
    for i=1:length(ns5SeqChannels)
        currVisCh = ns5VisChannels(i); % all filenames are based on visual channel numbers
        saveFilename = ['lfp' num2str(currVisCh) '_dec.mat'];
        saveFullPath = [processedDataDir saveFilename];
        if exist(saveFullPath, 'file') == 0
            % File Doesn't Exist. Add to "channels" list
            channelIndToProcess = [channelIndToProcess i]; 
        end
    end
end

%% NOW load one "original" 30000Hz LFP at a time and decimate it by a factor of 15
%  Note that we only care about visual channel numbers now - as we are no longer interacting with the openNSx function
for currVisCh=ns5VisChannels(channelIndToProcess) % XXX don't use "ch" as the loaded file contains the variable "ch"
    loadFilename = ['lfp' num2str(currVisCh) '_original.mat'];
    loadFullPath = [processedDataDir loadFilename];
    load(loadFullPath); % this will introduce lfp, Fs, metaTags, numPackets, visChNum, seqChNum
    
    % XXX insert here:
    % save the file in NS5 format somehow
    % filter the data using wavefilter, and resave the data in NS5 format -
    % this can then be used to process the file in plexon later
    
    % Decimation - NOW ONLY USING DECIMATION METHOD 1
    % Try 3 different ways of decimating by a factor of 15
    % 1. Using the matlab routine decimate - this uses chebyshev by default
    %    This is an IIR, that makes use of filtfilt so should be a 0 phase filt
    % 2. By hand - call my filterLFP function passing in 0 and 1000, then pick
    %    points 15:15:end of the data
    % 3. Don't low pass filter at all, just pick points 10:10:end

    %% TIME AXES
    FsFull = Fs;
    timeAxisFull = [1:length(lfp)] / FsFull;
    FsDec = FsFull / decFactor;
    %% DECIMATION STRATEGY 1 - matlab routine, it will lowpass filter below
    %% 0.8*Nyquist - so 800 Hz when converting to a sampling rate of 2000 Hz
    %lfpDec1 = decimate(double(lfp), decFactor);
    %lfpDec1 = resample(double(lfp(:)),1, decFactor);
    % timeAxisDec1 is special, and the following is based on observation of the
    % decimate.m file
    nd = length(lfp);
    r = decFactor;
    nout = ceil(nd/r);
    nbeg = r - (r*nout - nd);
    timeAxis = timeAxisFull(nbeg:r:nd);
    clear timeAxisFull
    
    %% DECIMATION STRATEGY 2 - butterworth filtfilt + pick every Nth pt
    %lfpFiltFilt = filterLFP(lfpData.Data, FsFull, 0.01, 800, 2, 0);
    %lfpDec2 = lfpFiltFilt(decFactor:decFactor:length(lfpData.Data));
    
    %% DECIMATION STRATEGY 3 - No low pass filter + pick every Nth pt
    %lfpDec1 = lfpData.Data(decFactor:decFactor:length(lfpData.Data));
    %lfpDec3 = lfpData.Data(decFactor:decFactor:length(lfpData.Data));

    Fs  = FsDec;
    lfp = round(lfpDec1); % round to the nearest integer to reduce file size
    clear lfpDec1
    % Calcualate the LFP stats necessary to compute the composite standard deviation
    % See http://www.burtonsys.com/climate/composite_standard_deviations.html
    clear lfpStats;
    lfpStats.mean = nanmean(lfp);
    lfpStats.var = nanvar(double(lfp));
    lfpStats.std = sqrt(lfpStats.var);
    lfpStats.N = length(lfp);
    lfpStats.nonNanN = length(find(isnan(lfp) == 0));
    lfpStats.errSumSq = lfpStats.var * (lfpStats.N - 1); % this is important for calculating the overall standard deviation
    lfpStats.rms = sqrt(nanmean(lfp.^2));
    lfpStats.max = nanmax(lfp);
    lfpStats.min = nanmin(lfp);
        
    saveFilename = ['lfp' num2str(currVisCh) '_dec.mat'];
    saveFullPath = [processedDataDir saveFilename];
    disp(['Saving Decimated LFP ' num2str(currVisCh)]); tic;
    numPackets = length(lfp);
    save(saveFullPath, 'lfp', 'Fs', 'metaTags', 'numPackets', 'visChNum', 'seqChNum', 'timeAxis', 'lfpStats'); % the decimated file should contain the same info, PLUS timeAxis (timeAxis not needed here as we are not skipping anything)
    toc;
    clear lfp lfpDec1 timeAxis Fs FsDec metaTags numPackets visChNum seqChNum lfpStats
end

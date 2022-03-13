%% GOALS of this script:
% - read in session list from session_list_human_lfp_spike_relationship ... DONE
% - convert NS5 to .mat DONE
% - stores metadata about each NS5 channel's signal stats: mean/variance/etc ... DONE
% - reads the corresponding NEVs and assigns a ms-PRECISE start/end time ... NOT POSSIBLE
% - filters and calls processLFPMinMax for each of the filter ranges needed
% - calculates the power spectrum over 2 second bins. Stores the signal
%   stats for each frequency (mean/variance/etc so that we can reconstruct
%   the overall mean/variance/etc for overall data).
% - plots the power spectrum for easy reference.
% - converts associated EDF file to .mat file, and stores start/end/mean/var
% - TRIGGER alignment: figures out the NS5 to EDF mapping and stores it.
% - [FUTURE: Extracts individual wave information for each channel]
% - So now we have all the data we need ready. The only thing left to do is
%   to analyze stretches of data and combine it with spikes.
% - Anything that is more seizure centric?
% - Once the LFP is processed like this (and the spikes are sorted), I can
%   also process a task using a custom written script. For task data I will
%   have to recalculate the triggered averages. 

%% 1. Source the session list, define sessions to analyze
% session_list_human_lfp_spike_relationship_sfn_onwards
% session_list_human_seizures_FOR_PAPER
% session_list_human_eyes_open_eyes_closed
session_list_human_seizures_NEW_FOR_ELLEN_SORTED_DATA_2017_TJ
ses = [6];

%processOrConcat = 1; % 1 = Process data & convert to .mat only; 2 = Concatenate data only; 3 = Both process and concatenate

%Tibin for pt RIHE1
processOrConcat = 3; % 1 = Process data & convert to .mat only; 2 = Concatenate data only; 3 = Both process and concatenate
processECoG = 0; % Whether to also process the associated ECoG File (convert each channel to .mat and generate the ns5-edf map).
checkBeforeConverting = 1; % If set to 1, it will NOT overwrite preexisting processed/converted data

filterBands = [];% DON'T USE THIS. Use the loops central_preprocessing function instead  [0.5 1; 1 2; 2 4; 4 6; 6 8; 8 10; 10 12; 12 16; 16 20; 20 32; 32 56; 56 64; 64 96;];
% Ideal Filter List
%  0.5 -   1 % Really Low Delta
%    1 -   2 % Low Delta
%    2 -   4 % High Delta
%    4 -   6 % Really Low Theta
%    6 -   8 % Low Theta
%    8 -  10 % High Theta / Low Alpha
%   10 -  12 % Alpha
%   12 -  16 % High Alpha
%   16 -  20 % Low Beta
%   20 -  32 % Beta
%   32 -  56 % Low Gamma
%   56 -  64 % Cover 60 Hz noise range - Call this nothing but still plot it
%   64 -  96 % High Gamma
%   96 - 256 % Ripple?
%  256 - 512 % MUA?
%filterBands = [];

NS5_PROCESSED_DEFAULT_DIR = 'J:/processed_data'; %   % /subject/ns5filename/. All data relative to this

%dataDir='/home/tibintj/turboHome/FOR_TIBIN_otherPts/seizures/RIHE1/NS5';
dataDir='/home/tibintj/turboHome/FOR_TIBIN_otherPts/seizures';

if processOrConcat == 1 || processOrConcat == 3
    for currSes = ses
        %% ECOG PROCESSING: First handle the ECoG file if there is one
        if processECoG == 1
            if exist('pop_biosig') ~= 2
                eeglab
            end
            for currEDFInd = 1:dataInfo(currSes).numEDFFiles % LOOP AROUND ALL NS5 FILES
                currEDFFilename = dataInfo(currSes).edfFilenames{currEDFInd};
                if ~isempty(currEDFFilename)
                    humanConvertEDFtoMat(checkBeforeConverting, dataInfo(currSes).subject, currEDFFilename, dataInfo(currSes).edfChannels);
                end
            end
        end
%humanConvertEDFtoMat(checkBeforeConverting, subject, edfFilename, edfSeqChannels, originalDataDir, processedDataDir)
        %% NS5 PROCESSING: Make sure that the relevant NS5 files have been converted to, for example, lfp6_original.mat and lfp6_dec.mat
        for currNS5Ind = 1:dataInfo(currSes).numNS5Files % LOOP AROUND ALL NS5 FILES
            currNS5Filename = dataInfo(currSes).ns5Filenames{currNS5Ind};
            humanConvertNS5toMatSeizures(checkBeforeConverting, dataInfo(currSes).subject, currNS5Filename, ...
                dataInfo(currSes).ns5VisChannels, dataInfo(currSes).ns5SeqChannels,dataDir,dataDir);
            if ~isempty(dataInfo(currSes).ns5TrigVisChannel)
                humanConvertNS5toMatSeizures(checkBeforeConverting, dataInfo(currSes).subject, currNS5Filename, ...
                    dataInfo(currSes).ns5TrigVisChannel, dataInfo(currSes).ns5TrigSeqChannel,dataDir,dataDir);
            end
        end
        
        %% Now get the Mean and Std for the ENTIRE session (i.e. averaged across all the ns5Files used to make up this session
        %clear overallLFPStats % This is a cell array that stores the overall LFP stats for each channel
        %for currCh = dataInfo(currSes).ns5VisChannels % LOOP AROUND ALL CHANNELS (because the humanGetOverallLFPStats will take care of looping around all files
         %   useDecimatedFiles = 1; % whether to work on 30000 or 2000 Hz files
         %   overallLFPStats{currCh} = humanGetOverallLFPStats(dataInfo(currSes).subject, currCh, dataInfo(currSes).numNS5Files, dataInfo(currSes).ns5Filenames, useDecimatedFiles);
        %end
  
        %% TRIGGER ALIGNMENT IF USING ECOG
        if processECoG == 1
            if exist('pop_biosig') ~= 2
                eeglab
            end
            if dataInfo(currSes).numNS5Files ~= dataInfo(currSes).numEDFFiles
                error('Mistmatching number of NS5 and EDF files - cannot align triggers');
            end
            for currInd = 1:dataInfo(currSes).numNS5Files % LOOP AROUND ALL NS5 FILES
                currNS5Folder = dataInfo(currSes).ns5Filenames{currInd};
                currEDFFolder = dataInfo(currSes).edfFilenames{currInd};
                trigNS5Filename = ['lfp' num2str(dataInfo(currSes).ns5TrigVisChannel) '_dec.mat'];
                trigEDFFilename = ['ecog_' dataInfo(currSes).edfTrigVisChannel '_original.mat'];
                
                if isempty(currEDFFolder)
                    warning(['Cannot align with NS5 file ' currNS5Filename ' - no matching EDF file provided']);
                    continue;
                end
                alignEDFNS5TriggersFromMat(dataInfo(currSes).subject, currEDFFolder, currNS5Folder, trigEDFFilename, trigNS5Filename, 1);
            end
        end
        
        %% MAIN PRE-PROCESSING LOOPS
        % Now for EACH FILE, and EACH CHANNEL, PRE-PROCESS the data to generate filtered data
        %         for currNS5Ind = 1:dataInfo(currSes).numNS5Files % LOOP AROUND ALL NS5 FILES
        %             currNS5Filename = dataInfo(currSes).ns5Filenames{currNS5Ind};
        %
        %             for currCh = dataInfo(currSes).ns5VisChannels % LOOP AROUND ALL CHANNELS
        %                 % 1. Load the LFP Data
        %                 loadPrefix   = ['lfp' num2str(currCh) '_dec'];
        %                 loadFilename = [loadPrefix '.mat'];
        %                 currProcessedDataDir = sprintf('%s/%s/%s/', NS5_PROCESSED_DEFAULT_DIR, dataInfo(currSes).subject,  currNS5Filename);
        %                 loadFullPath = [currProcessedDataDir loadFilename];
        %                 loadFullPrefix = [currProcessedDataDir loadPrefix];
        %                 lfpData = load(loadFullPath); % these var are now available: 'lfp', 'Fs', 'timeAxis', 'metaData', 'numPackets', 'ch'
        %
        %                 for currFilterInd = 1:size(filterBands,1)
        %                     % 2. This is a preprocessing script. So filter the LFP data in all the
        %                     % required bands. Just like in omarConvertNS5toMat, save a metadata
        %                     % file that contains the mean, std, err, etc of the filtered data trace
        %                     % 3. Now call processLFPMinMax (both steps 2 and 3 can be done here).
        %                     %    Call this file humanProcessLFPMinMax and it will be responsible
        %                     %    for first checking to see if the filtered LFP file exists. Then it
        %                     %    will check if the meta data exists. The filtered LFP file needs to
        %                     %    be in Volts, not Zscore. Attach the suffix volts to make this
        %                     %    clear. Min-max detection should still run on the Zscored data,
        %                     %    because otherwise I wouldn't really know what to set the threshold
        %                     %    value to. Keep the threshold at 0.2. If I need to check raw
        %                     %    amplitude of the filtered signal I can always just load the
        %                     %    filtered and saved trace.
        %                     currLowFreq = filterBands(currFilterInd,1);
        %                     currHighFreq = filterBands(currFilterInd,2);
        %                     disp(['About to process: File ' currNS5Filename ', Ch ' num2str(currCh) ', Band ' num2str(currLowFreq) '-' num2str(currHighFreq) '.']);
        %                     ripCont = humanProcessLFPMinMax(lfpData.timeAxis, lfpData.lfp, lfpData.Fs, currLowFreq, currHighFreq, loadFullPrefix, 1); % Note the ,1 addition that tells it to save the filteredLFP
        %                 end % End loop around filter bands
        %             end % End loop around channels
        %         end % End loop around NS5 files that are being analyzed together
    end % End loop around ses
end

if processOrConcat == 2 || processOrConcat == 3
    for currSes = ses
	dataInfo(currSes).subject
        for i=dataInfo(currSes).ns5VisChannels
        currOutputFilename = sprintf('%s_ch%d', dataInfo(currSes).descriptiveFilename, i);
	disp('here')
        omarCombineLFPChanIntoNS5(dataInfo(currSes).subject, i, dataInfo(currSes).numNS5Files, dataInfo(currSes).ns5Filenames, currOutputFilename, dataInfo(currSes).descriptiveFilename, 0,dataDir);
        end
    end
end
% FOR propagation analysis.
% I care about what's happening on a spike by spike basis. For each spike what were the LFP xcorr lags
% 1. So for each spike, take the xcorr centered on the spike for +- 100 ms. Store the time and the offset for each electrode.
% 2. So I will end up with a 2d matrix: rows are spike numbers, columns are electrode numbers. And the data contains max-lag-offset.
% 3. For each electrode, I also know its relative x,y location compared to the "local" (spike) electrode.
% 4. Now pick one non-local electrode. Find all the cases in which it preceeded (in time) the spike, and all the cases in which it came after
% 5. Plot the spike-triggered LFP in each case.
% 6. Plot the phase histogram in each case.
% 7. Calculate Mean Phase + Kappa + Rayleigh Z

%% HOW TO PROCESS A CONCATENATED LFP FOR SINGLE CYCLE ANALYSIS
% Scenario A - you want to use an adjusted threshold for all the files
% 1. First filter all the concatenated NS5s in a one filter range using humanProcessLFPMinMax, passing onlyFilterDontMinMax = 1
% 2. Then call getOverallFilteredLFPStats to get the means and SDs for each file
% 3. For the OVERALL concat file, 1Z = (x - overallMean)/overallSD. Therefore for 0.2Z, x = (0.2*overallSD) + overallMean
% 4. So for a single file we want to filter at a single file value of x
% 5. This value of x corresponds to a single file Z value of thresh = ([(0.2*overallSD) + overallMean] - singleMean) / singleSD
% 6. The adjusted threshold for each file, peakDetectThresh is thus ([(0.2*overallSD) + overallMean] - singleMean) / singleSD
% 7. Then call humanProcessLFPMinMax again, this time to calculate the minmax (i.e. set onlyFilterDOntMinMax to 0) and pass in the adjusted peakDetectThreshZ
% 8. This does not solve the main problem though, which is that some files have much lower or higher overall gamma values and zscoring individually gives z-score cycle amp distributions that are not comparable
% 9. Since the whole point is to get the zscores the same, I should pass in the overallMean and overallSD into humanProcessMinMax
% 10. Then this controls which type of Zscoring (single vs overall is used). overallZ will be added to the filename for those files
% 11. Thus there is no post-processing to do in this function

%% DEFAULTS
% CONCATENATED_FILE_OUTPUT_DEFAULT_DIR = '../TEMP_ANALYSIS_DEPOT/concatenated_ns5s_for_clustering/'; % /outputParentDir/outputFilename .Where to store the data
% NS5_PROCESSED_DEFAULT_DIR = '../TEMP_ANALYSIS_DEPOT/processed_data'; %   % /subject/ns5filename/. All data related to this
% CELLS_PROCESSED_DEFAULT_DIR = '../TEMP_ANALYSIS_DEPOT/sorted_units/concatenated'; % /subject/
% FIGURES_DEFAULT_DIR = '../TEMP_ANALYSIS_DEPOT/figures/'; %   % /subject/ns5filename/. All data related to this
CONCATENATED_FILE_OUTPUT_DEFAULT_DIR = 'K:/concatenated_ns5s_for_clustering/'; % /outputParentDir/outputFilename .Where to store the data
NS5_PROCESSED_DEFAULT_DIR = 'K:/processed_data'; %   % /subject/ns5filename/. All data related to this
CELLS_PROCESSED_DEFAULT_DIR = 'H:/sorted_units/concatenated'; % /subject/   % leaving this on G: as the units are not that large but are accessed often, so best to keep it fast
FIGURES_DEFAULT_DIR = 'K:/figures/'; %   % /subject/ns5filename/. All data related to this

%% SETTINGS
session_list_human_lfp_spike_relationship_sfn_onwards;
% session_list_human_seizures_FOR_PAPER;
ses = [1]; % NOT DONE (PROBLEM CHANNELS): ses1(MG29long) ch81&82&87, ses3(MG49main) ch50 (anesthesia), ses9(MG67task) ch80
% BAD CHANNEL LIST:
% MG29(ses1): ch81,82,87
% MG29(ses2): ch1,2,3,4,5 (anes),
% MG49(ses3): ch50 (anes) ch89 (disconnected)
% MG67(ses9): ch80
% T1, badCh=34,36,38,40,42,70, notRec=ch6&44

lfpChToProcess = [1 43 51]; % LFP CHANNELS, NOT spikes
maxElecDist = [400]; % this will process all lfp channels 400 microns from the channels listed in lfpChToProcess
ecogChToProcess = [];
processSingleCycle = 1;
overwriteMinMax = 0;
calcSpec = 1;
plotAndSaveSpecFig = 0;

preFilterLFP = 1; % This will filter the LFP between 0.5 and 900 Hz before calculating the spectrogram
windowSize = 1; % NOTE: using 2 second window to allow 0.5 Hz to be better represented as well
tapers = [3 5]; % if using 0.5 sec windows, switch to [1 1] to improve freq resolution
fpass = [0.5 250];
numLogBins = 50;
lfpFreq = [0.5 4; 4 8; 8 12; 12 18; 18 25; 25 30; 30 55; 65 110]; %[0.5 4; 4 8; 8 12; 30 55; 40 80; 0.5 1; 1 2; 2 3; 3 4; 4 5; 5 6; 6 7; 7 8; 8 9; 9 10; 10 11; 11 12; 12 14; 14 16; 16 18; 18 20; 20 22; 22 24; 24 26; 26 28; 28 30; 30 34; 34 38; 38 42; 42 46; 46 50; 50 54; 54 58; 58 62; 62 66; 66 70;]% 70 74; 74 78; 78 82; 82 86; 86 90];%[0.5 4; 4 8; 8 12; 12 18; 18 25; 25 32; 32 55; 65 110];
lfpFreqName = {'Delta', 'Theta', 'Alpha', 'Low Beta', 'High Beta', 'Beta-Gamma-Boundary', 'Low Gamma', 'High Gamma'};
plotLFPInd = 1%[1 3 7];
plotLFPColorMap = [0.5 0.5 0.5; 0.8 0.2 0.2; 0.2 0.8 0.2];
numFreq = size(lfpFreq,1);

if processSingleCycle == 1
    %% LOOP AROUND SESSIONS
    for currSes = ses
        [arrayMap electrodeXY electrodeImp] = neuroportArrayData(dataInfo(currSes).subject);
        % If no LFP channels provided - then ignore distance and process all channels
        if isempty(lfpChToProcess) % if no distance, and no lfp channels provided, process them all.
            allLFPChToProcess = dataInfo(currSes).ns5VisChannels;
        else
            if isempty(maxElecDist) % if no elecDist provided
                allLFPChToProcess = lfpChToProcess;
            else
                allLFPChToProcess = [];
                for tempLFPCh = lfpChToProcess
                    currLFPChToProcess = neuroportArrayGetElecAtDist(tempLFPCh, 0, maxElecDist, arrayMap, electrodeXY);
                    allLFPChToProcess = [allLFPChToProcess currLFPChToProcess];
                end
                % now find the unique ones and the good ones
                allLFPChToProcess = intersect(unique(allLFPChToProcess), dataInfo(currSes).ns5VisChannels);
                disp(['LFPsToProcess=' num2str(lfpChToProcess) ', MaxDist=' num2str(maxElecDist) ', FINALLFPChToProcess=' num2str(allLFPChToProcess)]);
            end
        end
        
        %% LOOP AROUND ALL CHANNELS (to plot)
        for currCh = allLFPChToProcess
            
            %% LOOP AROUND FREQUENCIES
            % call this currMainFreqInd to indicate that it is a MAIN wrapper (because currFreq and currFreqInd are used in subroutines inside the loop around
            % NS5 files
            for currFreqInd = 1:numFreq
                currLowFreq  = lfpFreq(currFreqInd,1);
                currHighFreq = lfpFreq(currFreqInd,2);
                
                %% LOOP AROUND NS5 files
                % for each file call humanProcessLFPMinMax to process that lfp
                clear cycleInfo;
                for currNS5Ind = 1:dataInfo(currSes).numNS5Files
                    outputString = sprintf('PREPROCESSING FILTERED LFP: About to Process File %s, NS5Ind=%g, LFP=%g', dataInfo(currSes).ns5Filenames{currNS5Ind}, currNS5Ind, currCh); disp(outputString);
                    % Load the LFP
                    currProcessedDataDir = sprintf('%s/%s/%s/', NS5_PROCESSED_DEFAULT_DIR, dataInfo(currSes).subject, dataInfo(currSes).ns5Filenames{currNS5Ind});
                    loadFilename = ['lfp' num2str(currCh) '_dec.mat'];
                    loadFullPath = [currProcessedDataDir loadFilename];
                    lfpFilenamePrefix = [currProcessedDataDir 'lfp' num2str(currCh) '_dec'];
                    lfpData = load(loadFullPath);
                    humanProcessLFPMinMax(lfpData.timeAxis, lfpData.lfp, lfpData.Fs, currLowFreq, currHighFreq, lfpFilenamePrefix, 1, 1); % Note the ,1,1 addition that tells it to save the filteredLFP and onlyFilterDontMixMax (both are 1)
                    % ECoG Processing RIGHT HERE
                    % Load the mapNS5EDF file for this NS5 file. Then
                    % load/setup the ECoG data that corresponds to that
                    % file and only that part. Then call
                    % humanProcessLFPMinMax 
                end
                % Now all NS5s for this freq for this channel have been processed... so get the "overall" filtered stats
                % Then call the humanProcessMinMax function again - this time to calculate the minMax
                % Decide if the minMax should be stored as regular Zscore or
                % MZscore, or if I should just remove th
                overallLFPStats = humanGetOverallFilteredLFPStats(dataInfo(currSes).subject, currCh, dataInfo(currSes).numNS5Files, dataInfo(currSes).ns5Filenames, currLowFreq, currHighFreq, NS5_PROCESSED_DEFAULT_DIR);
                for currNS5Ind = 1:dataInfo(currSes).numNS5Files
                    outputString = sprintf('PREPROCESSING MIN MAX: About to Process File %s, NS5Ind=%g, LFP=%g', dataInfo(currSes).ns5Filenames{currNS5Ind}, currNS5Ind, currCh); disp(outputString);
                    % Load the LFP
                    currProcessedDataDir = sprintf('%s/%s/%s/', NS5_PROCESSED_DEFAULT_DIR, dataInfo(currSes).subject, dataInfo(currSes).ns5Filenames{currNS5Ind});
                    loadFilename = ['lfp' num2str(currCh) '_dec.mat'];
                    loadFullPath = [currProcessedDataDir loadFilename];
                    lfpFilenamePrefix = [currProcessedDataDir 'lfp' num2str(currCh) '_dec'];
                    lfpData = load(loadFullPath);
                    useOverallMeanStd = 1;
                    cycleInfo{currNS5Ind} = humanProcessLFPMinMax(lfpData.timeAxis, lfpData.lfp, lfpData.Fs, currLowFreq, currHighFreq, lfpFilenamePrefix, 1,0,overwriteMinMax, useOverallMeanStd, overallLFPStats.mean, overallLFPStats.std); % Note the ,1,0,1 addition that tells it to save the filteredLFP=1,onlyFilterDontMixMax=0,overwriteMinMax = 1
                end
            end
        end
    end
end

%% SEPARATE LOOP TO PROCESS ALL THE SPECTROGRAMS and SAVE SPECTRAL FIGURES
if calcSpec == 1
    for currSes = ses
        [arrayMap electrodeXY electrodeImp] = neuroportArrayData(dataInfo(currSes).subject);
        % If no LFP channels provided - then ignore distance and process all channels
        if isempty(lfpChToProcess) % if no distance, and no lfp channels provided, process them all.
            allLFPChToProcess = dataInfo(currSes).ns5VisChannels;
        else
            if isempty(maxElecDist) % if no elecDist provided
                allLFPChToProcess = lfpChToProcess;
            else
                allLFPChToProcess = [];
                for tempLFPCh = lfpChToProcess
                    currLFPChToProcess = neuroportArrayGetElecAtDist(tempLFPCh, 0, maxElecDist, arrayMap, electrodeXY);
                    allLFPChToProcess = [allLFPChToProcess currLFPChToProcess];
                end
                % now find the unique ones and the good ones
                allLFPChToProcess = intersect(unique(allLFPChToProcess), dataInfo(currSes).ns5VisChannels);
                disp(['LFPsToProcess=' num2str(lfpChToProcess) ', MaxDist=' num2str(maxElecDist) ', FINALLFPChToProcess=' num2str(allLFPChToProcess)]);
            end
        end
        
        %% LOOP AROUND ALL CHANNELS (to plot)
        for currCh = allLFPChToProcess
            %% LOOP AROUND NS5 files
            for currNS5Ind = 1:dataInfo(currSes).numNS5Files
                outputString = sprintf('PREPROCESSING SPECTRUM: About to Process File %s, NS5Ind=%g, LFP=%g', dataInfo(currSes).ns5Filenames{currNS5Ind}, currNS5Ind, currCh); disp(outputString);
                % Load the LFP
                currProcessedDataDir = sprintf('%s/%s/%s/', NS5_PROCESSED_DEFAULT_DIR, dataInfo(currSes).subject, dataInfo(currSes).ns5Filenames{currNS5Ind});
                loadFilename = ['lfp' num2str(currCh) '_dec.mat'];
                loadFullPath = [currProcessedDataDir loadFilename];
                lfpFilenamePrefix = [currProcessedDataDir 'lfp' num2str(currCh) '_dec'];
                lfpData = load(loadFullPath);
                % Calculate and Store Spectrogram
                specData = omarSpikeLFPRatePower(0, 0, 0, lfpData.timeAxis, lfpData.lfp, lfpData.Fs, lfpFilenamePrefix, windowSize, tapers, fpass, numLogBins, preFilterLFP);
                if plotAndSaveSpecFig == 0
                    continue;
                end
                
                %% CALCULATE POWER IN EACH FREQ BAND
                allFreqPower = NaN(numFreq,  length(specData.t));
                allFreqPowerZ = NaN(numFreq, length(specData.t));
                allFreqPowerMZ = NaN(numFreq, length(specData.t));
                for currFreq = 1:numFreq
                    currLowFreq = lfpFreq(currFreq,1);
                    currHighFreq = lfpFreq(currFreq,2);
                    currFreqInd = find(specData.f >= currLowFreq & specData.f <= currHighFreq);
                    currMeanPower = mean(specData.S(:,currFreqInd),2);
                    allFreqPower(currFreq, :) = currMeanPower;
                    allFreqPowerZ(currFreq, :) = zscoreLFP(currMeanPower);
                    allFreqPowerMZ(currFreq, :) = mzscore(currMeanPower); % Very linear correlation between median-based Zscore and mean-based Zscore. But median z-score is better
                    % because when there are extremely large z-score values
                end
                
                %% plot and save the spectrum in ../../figures/MG49/$descriptiveFilename$/lfp6_dataInfo(currSes).ns5Filenames{currNS5Ind}_spectrum
                SLog = log(specData.S);
                SLogSmooth = omarSmooth2D(SLog,5,1); % smoothxy: x=time, y=freq. So increase the first param to smooth more across time without blurring frequencies
                lw = 1.5;
                fig = figPrep([10 6], 15);
                %subFig = subplot(11,4,[1:3 5:7 9:11 13:15 17:19 21:23]);
                subFig = subplot(18,1,1:7);
                omarPcolor(specData.t/100, specData.f, (SLogSmooth'), subFig);
                shading flat;
                caxis([0 8])
                ylim([0.5 70]);
                xvals = get(gca, 'xlim');
                set(gca, 'xtick', [0:100/100:max(specData.t)]);
                set(gca, 'tickdir', 'out');
                set(gca, 'xaxislocation', 'top');
                % RAW LFP
                subplot(18,1,8:10);
                plot(lfpData.timeAxis/100, zscore(lfpData.lfp), 'linewidth', 0.5, 'color', 'k');
                set(gca, 'xtick', [0:100/100:max(specData.t)]);
                xlim(xvals);
                ylim([-9 9]);
                set(gca, 'xticklabel', {});
                ylabel('Median Z (MZ)')
                % LOWER BAND POWERS
                subplot(18,1,11:18)
                cmap = jet(numFreq);
                hold on
                for currPlotIndInd = 1:length(plotLFPInd)
                    currFreqInd = plotLFPInd(currPlotIndInd);
                    plot(specData.t/100, fastsmooth(allFreqPowerMZ(currFreqInd,:),5,3,1), 'linewidth', lw, 'color', plotLFPColorMap(currPlotIndInd, :));
                    legText{currPlotIndInd} = lfpFreqName{currFreqInd};
                end
                grid on;
                xlim(xvals);
                set(gca, 'xtick', [0:100/100:max(specData.t)]);
                ylim([-2 5])
                leg = legend(legText);
                set(leg, 'box', 'off')
                xlabel('Time (100 seconds)');
                ylabel('Median-based Zscore (MZ)');
                
                printDir = [FIGURES_DEFAULT_DIR dataInfo(currSes).subject '/' dataInfo(currSes).descriptiveFilename];
                if exist(printDir) ~= 7
                    mkdir(printDir)
                end
                printFilename = [printDir '/lfp' num2str(currCh), '_dec_' dataInfo(currSes).ns5Filenames{currNS5Ind} specData.spectrumLFPPostfix];
                print('-djpeg100', '-r200', [printFilename '.jpg']);
                print('-dpdf', '-r200', [printFilename '.pdf']);
                %%
                close(fig);
            end
        end
    end
end
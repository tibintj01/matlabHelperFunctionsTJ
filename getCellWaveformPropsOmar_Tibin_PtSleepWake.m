function [cellProp] = getCellWaveformPropsOmar_Tibin_PtSleepWake(subject,sessIdx,descriptiveFilename, ch, indWithinCh, nexLoadDir, ns5DirName,saveDir)
%computes features of good spike waveforms for given cell (addressed by subject, session description, channel number,
%index within channel, and data directory) and saves as mat file to saveDir (and optionally images of waveforms)

% nexLoadDir and saveDir are optional. If provided they should each be a full path
% nexLoadDir refers to the location of the nex file that will be read in to process the spikes
% since the "descriptive_filename" folder in processed_data/MG49 is so packed, store the spike properties in:
% 20110615-094530-023_024_025_026_030_031_032_plus_anesthesia_SPKPROP

disp(sprintf('Save directory: %s',saveDir))

tic

if(isstr(ch))
	ch=str2num(ch);
end
if(isstr(indWithinCh))
	indWithinCh=str2num(indWithinCh);
end
%if(isstr(szIdx))
%	szIdx=str2num(szIdx);
%end

%session_list_human_seizures_NEW_FOR_ELLEN_SORTED_DATA_2017_TJ
session_list_human_lfp_spike_relationship_for_students
%sessIdx=find(ismember({dataInfo(:).subject},subject));

%cellChannels=dataInfo(sessIdx).cellChannel;
%cellChIdxes=find(cellChannels==ch);
%cellChIdx=cellChIdxes(indWithinCh);

%szStartTime=dataInfo(sessIdx).szStartTime(szIdx);
%szEndTime=dataInfo(sessIdx).szEndTime(szIdx);

%if(szIdx>1)
%	prevSzEndTime=dataInfo(sessIdx).szEndTime(szIdx-1)
%else
%	prevSzEndTime=-900;
%end

%lowTimeLimit=prevSzEndTime+900;
%highTimeLimit=szStartTime-300;

%Tibin addition
cellID=sprintf('%d%s',ch,num2letter(indWithinCh));
%szID=subject;

plotFigures=0;
%plotFigures=1;

CONCATENATED_FILE_OUTPUT_DEFAULT_DIR = 'J:/concatenated_ns5s_for_clustering/'; % /outputParentDir/outputFilename .Where to store the data
NS5_PROCESSED_DEFAULT_DIR = 'J:/processed_data'; %   % /subject/ns5filename/. All data related to this
CELLS_PROCESSED_DEFAULT_DIR = '../../sorted_units/concatenated'; % /subject/   % leaving this on G: as the units are not that large but are accessed often, so best to keep it fast
FIGURES_DEFAULT_DIR = 'J:/figures/'; %   % /subject/ns5filename/. All data related to this

if exist('saveDir') == 0 % if no saveDir variable passed in
    saveDir = [NS5_PROCESSED_DEFAULT_DIR '/' subject '/' descriptiveFilename '_SPKPROP' '/']; %STORE HERE
end
currSuffix = convertNumberToLetter(indWithinCh);
%cellPropSaveFile = sprintf('%s%g%s_cell_properties.mat', saveDir, ch, currSuffix);
%cellPropSaveFile = fullfile(saveDir,sprintf('%g%s_cell_properties_%s_seizure%d.mat',ch, currSuffix,subject,szIdx));
cellPropSaveFile = fullfile(saveDir,sprintf('%g%s_cell_properties_%s.mat',ch, currSuffix,subject));

%% RETURN CELL PROPERTIES INFO if cell has already been converted
%if exist(cellPropSaveFile, 'file') == 2
%    % fill up cellProp struct and return it
%    cellProp = load([cellPropSaveFile]);
%    return;
%end

%% SETUP SAVE DIR
if exist(saveDir) ~= 7 % if the actual saveDir directory does not exist
    mkdir(saveDir)
end

%% FIGURE OUT WHAT FOLDER CONTAINS THE NEX FILE
if exist('nexLoadDir') == 0
    nexLoadDir = [CONCATENATED_FILE_OUTPUT_DEFAULT_DIR '/' subject '/' descriptiveFilename '/']; % LOAD FROM HERE
    nexSecondLoadDir = [CELLS_PROCESSED_DEFAULT_DIR '/' subject '/' descriptiveFilename '/' 'chan' num2str(ch) '/'];

    spikeFile = [nexLoadDir '/' descriptiveFilename '_ch' num2str(ch) '.nex'];
    if ~exist(spikeFile, 'file')
        spikeFile = [nexSecondLoadDir '/' descriptiveFilename '_ch' num2str(ch) '.nex'];
        disp('USING spikeFile from the sorted_units directory');
    else
        % XX NOT IMPLEMENTED YET
        disp('USING spikeFile from the concantenated_file_output_default directory');
    end

else % if a nexLoadDir has been passed into this function, just use that folder
    %spikeFile = [nexLoadDir '/' descriptiveFilename '_ch' num2str(ch) '.nex'];
    %fileInfo=dir([nexLoadDir '/' descriptiveFilename '_ch' num2str(ch) '-*.nex']);
   fileInfo=dir([nexLoadDir '/' descriptiveFilename '_ch' num2str(ch) '.nex']);
	if(strcmp(subject,'RIHE1'))
   		fileInfo=dir([nexLoadDir '/RIHE1_' descriptiveFilename '_ch' num2str(ch) '-01.nex']);
	end 
    spikeFile=fullfile(nexLoadDir,fileInfo(1).name)
    %spikeFile = [nexLoadDir '/' descriptiveFilename '_ch' num2str(ch) '-01.nex'];
end


%% LOAD THE SPIKES FROM NEX FILE
info = ft_read_spike(spikeFile);
Fs = info.hdr.FileHeader.Frequency;
chanLabel = ['Channel01' currSuffix '_wf'];
cellName = sprintf('cell%g%s', ch, currSuffix);
chanToLoad = chan(info, chanLabel);
spikeBins = double(info.timestamp{chanToLoad}); % XXX double is very important!!
spikeTimes = spikeBins / Fs;

%inTimeSpikeIdxes=find(spikeTimes>lowTimeLimit  & spikeTimes < highTimeLimit);
%spikeTimes=spikeTimes(inTimeSpikeIdxes);

numSpikes = length(spikeTimes);
waveforms = info.waveform{chanToLoad};

%fds


%waveforms=waveforms(:,inTimeSpikeIdxes);


%Tibin addition 11/15/17; load raw spikes from nex file
%ns5DirName=fullfile(nexLoadDir,'..','NS5');
rawSamplesBefore=9;
rawSamplesAfter=38;
%[rawWaveformsPlusMinus2ms,rawWaveforms]=getRawWaveformsOfCell(ns5DirName,descriptiveFilename,ch,spikeTimes,rawSamplesBefore,rawSamplesAfter);

%[rawWaveforms]=getRawWaveformsOfCell(ns5DirName,descriptiveFilename,ch,spikeTimes,rawSamplesBefore,rawSamplesAfter);
%waveforms = rawWaveforms;
descriptiveFilenameLFP=descriptiveFilename;
if(strcmp('20090216-184530-025_026_027_plus_anesthesia',descriptiveFilename))
        descriptiveFilenameLFP='20090216-184530-025_026_027_REMOVED_anesthesia';
end

if(strcmp(subject,'RIHE1'))
	descriptiveFilenameLFP=['RIHE1_' descriptiveFilename];
end


[filtWaveforms]=getFiltWaveformsOfCell(ns5DirName,descriptiveFilenameLFP,ch,spikeTimes,rawSamplesBefore,rawSamplesAfter);
waveforms = filtWaveforms;

%truncate num spikes to what could be retrieved from corresponding raw data
numSpikes=size(waveforms,2);
spikeTimes=spikeTimes(1:numSpikes);


%[filtWaveforms]=getWaveformsFiltAll(ns5DirName,descriptiveFilename,ch,spikeTimes,rawSamplesBefore,rawSamplesAfter);
%waveforms = filtWaveforms;


%figure; plot(info.waveform{chanToLoad})
%saveas(gcf,'testNEXwaves.tif')

%figure; plot(rawWaveforms)
%saveas(gcf,'testRawwaves.tif')


%% THIS WAS BEING DONE IN THE CODE FOR FIGURE 1 EPILEPSY PAPER to generate spikes



%% CALCULATIONS
% calculate spikeWidth and spikeAmp for each spike
disp(['processCellProperties: aboutToCalcSpikeWidth, N = ' num2str(numSpikes)]);
%tic;
dualSpikes = [];
singleSpikes = [];
dualP50Spikes = [];
singleP50Spikes = [];
dualP25Spikes = [];
singleP25Spikes = [];
dualP75Spikes = [];
singleP75Spikes = [];

% all of the following assumes a spike that starts with a negative-going phase.
spikeMinAmp = zeros(1,numSpikes); % this is just the MINIMA of the spike
spikeMinIndex = zeros(1,numSpikes); % the bin number where the spike MINIMA occurs. This can be used to align the spike to the minima
spikeMaxAmp = zeros(1,numSpikes); % the is the FIRST MAXIMA of the spike
spikeMaxIndex = zeros(1,numSpikes); % the bin number where the FIRST spike MAXIMA occurs. This can be used to align the spike to the FIRST MAXIMA

spikeT50Width = zeros(1,numSpikes); % this is the full width at half max of the trough (the width of the neg-going part)
spikeT75Width = zeros(1,numSpikes); % this is the full width at half max of the trough (the width of the neg-going part)
spikeT25Width = zeros(1,numSpikes); % this is the full width at half max of the trough (the width of the neg-going part)

spikeT35Width = zeros(1,numSpikes); % this is the full width at half max of the trough (the width of the neg-going part)
spikeT30Width = zeros(1,numSpikes); % this is the full width at half max of the trough (the width of the neg-going part)
spikeT20Width = zeros(1,numSpikes); % this is the full width at half max of the trough (the width of the neg-going part)
spikeT15Width = zeros(1,numSpikes); % this is the full width at half max of the trough (the width of the neg-going part)
spikeT10Width = zeros(1,numSpikes); % this is the full width at half max of the trough (the width of the neg-going part)

spikeT05Width = zeros(1,numSpikes); % this is the full width at half max of the trough (the width of the neg-going part)
spikeT01Width = zeros(1,numSpikes); % this is the full width at half max of the trough (the width of the neg-going part)
spikeTMinus05Width = zeros(1,numSpikes); % this is the full width at half max of the trough (the width of the neg-going part)
spikeTMinus10Width = zeros(1,numSpikes); % this is the full width at half max of the trough (the width of the neg-going part)
spikeTMinus15Width = zeros(1,numSpikes); % this is the full width at half max of the trough (the width of the neg-going part)
spikeTMinus20Width = zeros(1,numSpikes); % this is the full width at half max of the trough (the width of the neg-going part)

% The next three are widths of the positive going phase - these are calculated using OVERALL MAX AMP. Later on I may also calculate them using the FIRST MAX AMP. 
% We shall see if that is necessary
spikeP50Width = zeros(1,numSpikes); % this is the full width at half max of the peak (the width of the pos-going part)
spikeP25Width = zeros(1,numSpikes); % this is the full width at 25% of the max of the peak (the width of the pos-going part). This is good for distinguishing FS/RS cells
spikeP75Width = zeros(1,numSpikes); % 75% is closer to the peak than 25%. this is the full width at 75% of the max of the peak (the width of the pos-going part). This is good for distinguishing FS/RS cells

spikeTPAmp = zeros(1,numSpikes); % trough-to-first-peak amplitude
spikeTPWidth = zeros(1,numSpikes); % trough-to-first-peak duration
spikeTP1090Width = zeros(1,numSpikes); % trough-to-first-peak duration 10-90
spikeTP2080Width = zeros(1,numSpikes); % trough-to-first-peak duration 20-80

spikeOverallMaxAmp = zeros(1,numSpikes); % this is the OVERALL maxima of the spike (not just the first). This is relevant if there were multiple overlapping spikes.
spikeOverallMaxIndex = zeros(1,numSpikes); % the bin of the OVERALL maxima
spikeOverallMaxTPAmp = zeros(1,numSpikes); % the spike amp calculate from first trough to OVERALL peak
spikeOverallMaxTPWidth = zeros(1,numSpikes); % the spike amp calculated from first trough to OVERALL PEAK

% Not sure what this is good for.
spikeStartMinIndex = zeros(1,numSpikes); 
spikeStartMin = zeros(1,numSpikes);
spikeSPAmp = zeros(1,numSpikes);
spikeSPWidth = zeros(1,numSpikes);

% Integrated amp
spikeIntraAmp = zeros(1,numSpikes); % how is this calculated? TP or just max?
spikeIntraWidth = zeros(1,numSpikes); % how is this cal
spikeIntraMinBin=zeros(1,numSpikes); 

numBins = 48; % 48 time points = 1.6 ms
for spikeNum = 1:numSpikes
    if mod(spikeNum, 1000) == 0
       disp(['Processed ' num2str(spikeNum) ' spikes.']);
    end
    spikeLow = waveforms(:,spikeNum);

    % now spline the spike
    splineDt = .1; % XXX in bins (need to adjust based on sampling rate?)
    spikeHigh = spline(1:numBins, spikeLow, 1:splineDt:numBins);
    [spikeMinAmp(spikeNum) spikeMinIndex(spikeNum)] = min(spikeHigh); % should I restrict to first N bins? this should be a NEG value
    halfMinAmp = spikeMinAmp(spikeNum)/2;
    %Tibin addition
    t25Amp = spikeMinAmp(spikeNum)*0.25;
    t75Amp = spikeMinAmp(spikeNum)*0.75;
    
    aboveHalfIndices = find(spikeHigh <= halfMinAmp);

    %Tibin addition
    aboveT25Indices = find(spikeHigh <= t25Amp);
    aboveTMinus20Indices = find(spikeHigh <= spikeMinAmp(spikeNum)*(-0.20));
    aboveTMinus15Indices = find(spikeHigh <= spikeMinAmp(spikeNum)*(-0.15));
    aboveTMinus10Indices = find(spikeHigh <= spikeMinAmp(spikeNum)*(-0.10));
    aboveTMinus05Indices = find(spikeHigh <= spikeMinAmp(spikeNum)*(-0.05));
    aboveT01Indices = find(spikeHigh <= spikeMinAmp(spikeNum)*0.01);
    aboveT05Indices = find(spikeHigh <= spikeMinAmp(spikeNum)*0.05);
    aboveT10Indices = find(spikeHigh <= spikeMinAmp(spikeNum)*0.10);
    aboveT15Indices = find(spikeHigh <= spikeMinAmp(spikeNum)*0.15);
    aboveT20Indices = find(spikeHigh <= spikeMinAmp(spikeNum)*0.20);
    aboveT30Indices = find(spikeHigh <= spikeMinAmp(spikeNum)*0.30);
    aboveT35Indices = find(spikeHigh <= spikeMinAmp(spikeNum)*0.35);
    aboveT75Indices = find(spikeHigh <= t75Amp);

	aboveT75Indices=keepConsecIndices(aboveT75Indices);
	aboveT25Indices=keepConsecIndices(aboveT25Indices);
	aboveT35Indices=keepConsecIndices(aboveT35Indices);
	aboveT30Indices=keepConsecIndices(aboveT30Indices);
	aboveT20Indices=keepConsecIndices(aboveT20Indices);
	aboveT15Indices=keepConsecIndices(aboveT15Indices);
	aboveT10Indices=keepConsecIndices(aboveT10Indices);
	aboveT05Indices=keepConsecIndices(aboveT05Indices);
	aboveT01Indices=keepConsecIndices(aboveT01Indices);
	aboveTMinus05Indices=keepConsecIndices(aboveTMinus05Indices);
	aboveTMinus10Indices=keepConsecIndices(aboveTMinus10Indices);
	aboveTMinus15Indices=keepConsecIndices(aboveTMinus15Indices);
	aboveTMinus20Indices=keepConsecIndices(aboveTMinus20Indices);


    if (length(find(diff(aboveHalfIndices) > 1)) ~= 0)
        %warning(['multiple spikes in spike # ' num2str(spikeNum)]) ; % this is not working great - as many spikes seem to be going below the halfMin value twice during the post-rebound phase
        dualSpikes = [dualSpikes; spikeNum];
    else
        singleSpikes = [singleSpikes; spikeNum];
    end

	%Tibin addition
	aboveHalfIndices=keepConsecIndices(aboveHalfIndices);

    spikeT50Width(spikeNum) = (length(aboveHalfIndices) - 1) * splineDt; % in samples!!
    
    spikeT35Width(spikeNum) = (length(aboveT35Indices) - 1) * splineDt; % in samples!!
    spikeT30Width(spikeNum) = (length(aboveT30Indices) - 1) * splineDt; % in samples!!
    spikeT25Width(spikeNum) = (length(aboveT25Indices) - 1) * splineDt; % in samples!!
    spikeT20Width(spikeNum) = (length(aboveT20Indices) - 1) * splineDt; % in samples!!
    spikeT15Width(spikeNum) = (length(aboveT15Indices) - 1) * splineDt; % in samples!!
    spikeT10Width(spikeNum) = (length(aboveT10Indices) - 1) * splineDt; % in samples!!
    spikeT05Width(spikeNum) = (length(aboveT05Indices) - 1) * splineDt; % in samples!!
    spikeT75Width(spikeNum) = (length(aboveT75Indices) - 1) * splineDt; % in samples!!
    spikeT01Width(spikeNum) = (length(aboveT01Indices) - 1) * splineDt; % in samples!!
    spikeTMinus05Width(spikeNum) = (length(aboveTMinus05Indices) - 1) * splineDt; % in samples!!
    spikeTMinus10Width(spikeNum) = (length(aboveTMinus10Indices) - 1) * splineDt; % in samples!!
    spikeTMinus15Width(spikeNum) = (length(aboveTMinus15Indices) - 1) * splineDt; % in samples!!
    spikeTMinus20Width(spikeNum) = (length(aboveTMinus20Indices) - 1) * splineDt; % in samples!!


%% TROUGH2PEAK, first PEAK
    spikeHighSmooth = fastsmooth(spikeHigh, 20, 3, 1); % NEW OMAR 2014 addition: Replace with causalSmooth?
    diffSpikeHigh = [0 diff(spikeHighSmooth)];
    % check for the first point after the trough where the derivative becomes
    % positive again. this is the location of the first peak after the trough
    % XXX OMAR change May 17, 2011:
    % the next line was finding diffSpikeHigh values that were less than
    % 0 because the smooth and the spikeHigh are offset. Find the MinIndex
    % of the smooth and use this to search diffSpikeHigh instead
    [ignoreMe spikeMinIndexSmooth] = min(spikeHighSmooth);
    firstDecreasingIndex = find(diffSpikeHigh(spikeMinIndexSmooth+1:end) <= 0, 1);
    if isempty(firstDecreasingIndex) == 0
        spikeMaxIndex(spikeNum) = firstDecreasingIndex + spikeMinIndexSmooth - 1; % XXX errorcheck this
        if (spikeMaxIndex(spikeNum) - spikeMinIndex(spikeNum)) > 10 % at least 1 original bin long
            spikeMaxAmp(spikeNum) = spikeHigh(spikeMaxIndex(spikeNum));
            spikeTPAmp(spikeNum) = spikeMaxAmp(spikeNum) - spikeMinAmp(spikeNum);
            spikeTPWidth(spikeNum) = (spikeMaxIndex(spikeNum)-spikeMinIndex(spikeNum)+1) * splineDt;
            % also store the 10-90% rise times
            amp10 = spikeMinAmp(spikeNum)+(spikeMaxAmp(spikeNum)-spikeMinAmp(spikeNum))*0.1;
            amp90 = spikeMinAmp(spikeNum)+(spikeMaxAmp(spikeNum)-spikeMinAmp(spikeNum))*0.9;
            if spikeNum == 13
                here = 1
            end
            [ampBins] = interp1_same_x_allowed(spikeHigh(spikeMinIndex(spikeNum):spikeMaxIndex(spikeNum)), spikeMinIndex(spikeNum):spikeMaxIndex(spikeNum), [amp10 amp90], 'nearest');
            spikeTP1090Width(spikeNum) = (ampBins(2)-ampBins(1)+1) * splineDt;
            % also store the 20-80% rise times
            amp20 = spikeMinAmp(spikeNum)+(spikeMaxAmp(spikeNum)-spikeMinAmp(spikeNum))*0.2;
            amp80 = spikeMinAmp(spikeNum)+(spikeMaxAmp(spikeNum)-spikeMinAmp(spikeNum))*0.8;
            [ampBins] = interp1_same_x_allowed(spikeHigh(spikeMinIndex(spikeNum):spikeMaxIndex(spikeNum)), spikeMinIndex(spikeNum):spikeMaxIndex(spikeNum), [amp20 amp80], 'nearest');
            spikeTP2080Width(spikeNum) = (ampBins(2)-ampBins(1)+1) * splineDt;
        else
            % what to do if no minima is found?!
            % XXX keep a list of these NaN indices
            spikeMaxIndex(spikeNum) = NaN;%length(spikeHigh);
            spikeMaxAmp(spikeNum) = NaN;
            spikeTPAmp(spikeNum) = NaN;
            spikeTPWidth(spikeNum) = NaN;
            spikeTP1090Width(spikeNum) = NaN;
            spikeTP2080Width(spikeNum) = NaN;
        end
    else
        % what to do if no minima is found?!
        % XXX keep a list of these NaN indices
        spikeMaxIndex(spikeNum) = NaN;%length(spikeHigh);
        spikeMaxAmp(spikeNum) = NaN;
        spikeTPAmp(spikeNum) = NaN;
        spikeTPWidth(spikeNum) = NaN;
        spikeTP1090Width(spikeNum) = NaN;
        spikeTP2080Width(spikeNum) = NaN;
    end
%% TROUGH2PEAK, overall PEAK
    if spikeMinIndex(spikeNum) == length(spikeHigh) % make sure the trough was not the last point in the waveform
        spikeOverallMaxIndex(spikeNum) = NaN;
        spikeOverallMaxAmp(spikeNum) = NaN;
        spikeOverallMaxTPAmp(spikeNum) = NaN;
        spikeOverallMaxTPWidth(spikeNum) = NaN;
    else
        [spikeOverallMaxAmp(spikeNum) spikeOverallMaxIndex(spikeNum)] = max(spikeHigh(spikeMinIndex(spikeNum)+1:end));
        spikeOverallMaxIndex(spikeNum) = spikeOverallMaxIndex(spikeNum) + spikeMinIndex(spikeNum);
        spikeOverallMaxTPAmp(spikeNum) = spikeOverallMaxAmp(spikeNum) - spikeMinAmp(spikeNum);
        spikeOverallMaxTPWidth(spikeNum) = (spikeOverallMaxIndex(spikeNum)-spikeMinIndex(spikeNum)+1) * splineDt;
    end
%% CALCULATE the width of the PEAK (of the positive going phase)
    %  Use the overallMaxAmp to calculate the 50%, 25%, 75% max values and calculate the durations at those values.
    %  I think 75% will be best for FS/RS differentiation
    %  Start with 50%
    halfOverallMaxAmp = spikeOverallMaxAmp(spikeNum)*0.5;
    aboveHalfIndices = find(spikeHigh >= halfOverallMaxAmp);
    if (length(find(diff(aboveHalfIndices) > 1)) ~= 0)
        %warning(['multiple P50 peaks in spike # ' num2str(spikeNum)]) ;
        dualP50Spikes = [dualP50Spikes; spikeNum];
    else
        singleP50Spikes = [singleP50Spikes; spikeNum];
    end
    spikeP50Width(spikeNum) = (length(aboveHalfIndices) - 1) * splineDt; % in samples!!
    %  25%
    overallMaxAmp25 = spikeOverallMaxAmp(spikeNum)*0.25;
    aboveHalfIndices = find(spikeHigh >= overallMaxAmp25);
    if (length(find(diff(aboveHalfIndices) > 1)) ~= 0)
        %warning(['multiple P25 peaks in spike # ' num2str(spikeNum)]) ;
        dualP25Spikes = [dualP25Spikes; spikeNum];
    else
        singleP25Spikes = [singleP25Spikes; spikeNum];
    end
    spikeP25Width(spikeNum) = (length(aboveHalfIndices) - 1) * splineDt; % in samples!!
    %  75%
    overallMaxAmp75 = spikeOverallMaxAmp(spikeNum)*0.75;
    aboveHalfIndices = find(spikeHigh >= overallMaxAmp75);
    if (length(find(diff(aboveHalfIndices) > 1)) ~= 0)
        %warning(['multiple P75 peaks in spike # ' num2str(spikeNum)]) ;
        dualP75Spikes = [dualP75Spikes; spikeNum];
    else
        singleP75Spikes = [singleP75Spikes; spikeNum];
    end
    spikeP75Width(spikeNum) = (length(aboveHalfIndices) - 1) * splineDt; % in samples!!


%% START2PEAK, compensate for non-zero starting values of spikes - skipping for now (May 16, 2014)
% if the peak is at value 1, this is a useless statistic, so set to NaN.
% otherwise, go ahead and calculate it.
%     if spikeMaxIndex(spikeNum) == 1
%         spikeStartMinIndex(spikeNum) = NaN;
%         spikeStartMin(spikeNum) = NaN;
%         spikeSPAmp(spikeNum) = NaN;
%         spikeSPWidth(spikeNum) = NaN;
%     else
%         [spikeStartMin(spikeNum) spikeStartMinIndex(spikeNum)] = min(spikeHigh(1:spikeMaxIndex(spikeNum)));
%         spikeSPAmp(spikeNum) = spikeAmp(spikeNum) - spikeStartMin(spikeNum);
%         halfMaxSPAmp = spikeSPAmp(spikeNum)/2;
%         aboveHalfSPIndices = find(spikeHigh >= halfMaxSPAmp);
%         spikeSPWidth(spikeNum) = (length(aboveHalfSPIndices) - 1) * splineDt; % in samples!!
%     end

%% INTEGRATED SPIKE (CUMSUM), width at half max
    intraSpikeHigh = cumsum(spikeHigh) * -1; % * -1 added May 16, 2014. integrate and multiply by -1 to make positive going
    spikeIntraAmp(spikeNum) = max(intraSpikeHigh);
    halfMaxAmp = spikeIntraAmp(spikeNum)/2;
    aboveHalfIndices = find(intraSpikeHigh >= halfMaxAmp);
    spikeIntraWidth(spikeNum) = (length(aboveHalfIndices) - 1) * splineDt; % in samples!!
   
    %Tibin addition 11/1/17
    [dummyMin, spikeIntraMinBin(spikeNum)]=min(intraSpikeHigh);

	%Tibin addition Nov.4,2017

	%residualDataFileName=sprintf('/home/tibintj/avgResidualFeatures%s_%s.mat',cellID,szID);
	%if(exist(residualDataFileName,'file'))
	%	residualData=load(residualDataFileName);
	%	spikeAvgTroughRiseDiff=-residualData.avgOvershootTrough;
	%	spikeAvgPeakDecayDiff=residualData.avgOvershoot;
	%else
		%fds
	%	spikeAvgTroughRiseDiff=NaN;
	%	spikeAvgPeakDecayDiff=NaN;
	%end

    %residualFileName=sprintf('avgResidualFeatures%s_%s.mat',cellID,subject);
    %minResidueAfterPeak=
 
    %     if spikePTWidth(spikeNum) < 3
    %         spikePChip = interp1(1:48,spikeLow,1:splineDt:48,'pchip');
    %         figure; hold on;
    %         plot([1:10:480], spikeLow, 'b');
    %         plot(spikeHigh, 'k');
    %         plot(spikePChip, 'r');
    %     end
end
%toc;
% IMPORTANT - convert spikeMaxIndex into low resolution (1-48) bin counts.
% Currently it is from the spikeHigh (bins 1-471)
spikeMinIndex = (spikeMinIndex - 1/splineDt + 1) .* splineDt;
spikeMaxIndex = (spikeMaxIndex - 1/splineDt + 1) .* splineDt;
spikeOverallMaxIndex = (spikeOverallMaxIndex - 1/splineDt + 1) .* splineDt;
%spikeStartMinIndex = (spikeStartMinIndex - 1/splineDt + 1) .* splineDt;

nonNanSpikes = find(isnan(spikeTPWidth) == 0); % XXX save this variable for future lookup

% get the subset of spikes with trough between 6 to 16 binswithin the singlespikes... these
% are the best single,well-alligned spikes and should be used to estimate
% the average wave-form and the meanWidth and meanAmp
tempGoodSpikes = find(spikeMinIndex(singleSpikes) >= 6 & spikeMinIndex(singleSpikes) <= 16);
goodSpikes = singleSpikes(tempGoodSpikes);
clear tempGoodSpikes
clear spikeOverallMax spikeStartMax % not saving these

%% MEANS
% min amp
spikeMinAmpStats.allMean    = nanmean(spikeMinAmp);
spikeMinAmpStats.singleMean = nanmean(spikeMinAmp(singleSpikes));
spikeMinAmpStats.goodMean   = nanmean(spikeMinAmp(goodSpikes));
spikeMinAmpStats.allStd     = nanstd(spikeMinAmp);
spikeMinAmpStats.singleStd  = nanstd(spikeMinAmp(singleSpikes));
spikeMinAmpStats.goodStd    = nanstd(spikeMinAmp(goodSpikes));
% max amp
spikeMaxAmpStats.allMean    = nanmean(spikeMaxAmp);
spikeMaxAmpStats.singleMean = nanmean(spikeMaxAmp(singleSpikes));
spikeMaxAmpStats.goodMean   = nanmean(spikeMaxAmp(goodSpikes));
spikeMaxAmpStats.allStd     = nanstd(spikeMaxAmp);
spikeMaxAmpStats.singleStd  = nanstd(spikeMaxAmp(singleSpikes));
spikeMaxAmpStats.goodStd    = nanstd(spikeMaxAmp(goodSpikes));
% half trough width
spikeT50WidthStats.allMean    = nanmean(spikeT50Width);
spikeT50WidthStats.singleMean = nanmean(spikeT50Width(singleSpikes));
spikeT50WidthStats.goodMean   = nanmean(spikeT50Width(goodSpikes));
spikeT50WidthStats.allStd     = nanstd(spikeT50Width);
spikeT50WidthStats.singleStd  = nanstd(spikeT50Width(singleSpikes));
spikeT50WidthStats.goodStd    = nanstd(spikeT50Width(goodSpikes));
% trough2peak amp
spikeTPAmpStats.allMean    = nanmean(spikeTPAmp);
spikeTPAmpStats.singleMean = nanmean(spikeTPAmp(singleSpikes));
spikeTPAmpStats.goodMean   = nanmean(spikeTPAmp(goodSpikes));
spikeTPAmpStats.allStd     = nanstd(spikeTPAmp);
spikeTPAmpStats.singleStd  = nanstd(spikeTPAmp(singleSpikes));
spikeTPAmpStats.goodStd    = nanstd(spikeTPAmp(goodSpikes));
% trough2peakh width
spikeTPWidthStats.allMean    = nanmean(spikeTPWidth);
spikeTPWidthStats.singleMean = nanmean(spikeTPWidth(singleSpikes));
spikeTPWidthStats.goodMean   = nanmean(spikeTPWidth(goodSpikes));
spikeTPWidthStats.allStd     = nanstd(spikeTPWidth);
spikeTPWidthStats.singleStd  = nanstd(spikeTPWidth(singleSpikes));
spikeTPWidthStats.goodStd    = nanstd(spikeTPWidth(goodSpikes));
% trough2peak10-90 width
spikeTP1090WidthStats.allMean    = nanmean(spikeTP1090Width);
spikeTP1090WidthStats.singleMean = nanmean(spikeTP1090Width(singleSpikes));
spikeTP1090WidthStats.goodMean   = nanmean(spikeTP1090Width(goodSpikes));
spikeTP1090WidthStats.allStd     = nanstd(spikeTP1090Width);
spikeTP1090WidthStats.singleStd  = nanstd(spikeTP1090Width(singleSpikes));
spikeTP1090WidthStats.goodStd    = nanstd(spikeTP1090Width(goodSpikes));
% trough2peak10-90 width
spikeTP2080WidthStats.allMean    = nanmean(spikeTP2080Width);
spikeTP2080WidthStats.singleMean = nanmean(spikeTP2080Width(singleSpikes));
spikeTP2080WidthStats.goodMean   = nanmean(spikeTP2080Width(goodSpikes));
spikeTP2080WidthStats.allStd     = nanstd(spikeTP2080Width);
spikeTP2080WidthStats.singleStd  = nanstd(spikeTP2080Width(singleSpikes));
spikeTP2080WidthStats.goodStd    = nanstd(spikeTP2080Width(goodSpikes));
% trough2peak overall min amp
spikeOverallMaxTPAmpStats.allMean    = nanmean(spikeOverallMaxTPAmp);
spikeOverallMaxTPAmpStats.singleMean = nanmean(spikeOverallMaxTPAmp(singleSpikes));
spikeOverallMaxTPAmpStats.goodMean   = nanmean(spikeOverallMaxTPAmp(goodSpikes));
spikeOverallMaxTPAmpStats.allStd     = nanstd(spikeOverallMaxTPAmp);
spikeOverallMaxTPAmpStats.singleStd  = nanstd(spikeOverallMaxTPAmp(singleSpikes));
spikeOverallMaxTPAmpStats.goodStd    = nanstd(spikeOverallMaxTPAmp(goodSpikes));
% peak2trough overall min width
spikeOverallMaxTPWidthStats.allMean    = nanmean(spikeOverallMaxTPWidth);
spikeOverallMaxTPWidthStats.singleMean = nanmean(spikeOverallMaxTPWidth(singleSpikes));
spikeOverallMaxTPWidthStats.goodMean   = nanmean(spikeOverallMaxTPWidth(goodSpikes));
spikeOverallMaxTPWidthStats.allStd     = nanstd(spikeOverallMaxTPWidth);
spikeOverallMaxTPWidthStats.singleStd  = nanstd(spikeOverallMaxTPWidth(singleSpikes));
spikeOverallMaxTPWidthStats.goodStd    = nanstd(spikeOverallMaxTPWidth(goodSpikes));
% spikeP50Width the width after 50% to the peak
spikeP50WidthStats.allMean = nanmean(spikeP50Width);
spikeP50WidthStats.singleMean = nanmean(spikeP50Width(singleSpikes));
spikeP50WidthStats.goodMean = nanmean(spikeP50Width(goodSpikes));
spikeP50WidthStats.singlePeakMean = nanmean(spikeP50Width(singleP50Spikes));
spikeP50WidthStats.allStd = nanstd(spikeP50Width);
spikeP50WidthStats.singleStd = nanstd(spikeP50Width(singleSpikes));
spikeP50WidthStats.goodStd = nanstd(spikeP50Width(goodSpikes));
spikeP50WidthStats.singlePeakStd = nanstd(spikeP50Width(singleP50Spikes));
% spikeP25Width the width after 25% to the peak
spikeP25WidthStats.allMean = nanmean(spikeP25Width);
spikeP25WidthStats.singleMean = nanmean(spikeP25Width(singleSpikes));
spikeP25WidthStats.goodMean = nanmean(spikeP25Width(goodSpikes));
spikeP25WidthStats.singlePeakMean = nanmean(spikeP25Width(singleP25Spikes));
spikeP25WidthStats.allStd = nanstd(spikeP25Width);
spikeP25WidthStats.singleStd = nanstd(spikeP25Width(singleSpikes));
spikeP25WidthStats.goodStd = nanstd(spikeP25Width(goodSpikes));
spikeP25WidthStats.singlePeakStd = nanstd(spikeP25Width(singleP50Spikes));
% spikeP75Width the width after 75% to the peak
spikeP75WidthStats.allMean = nanmean(spikeP75Width);
spikeP75WidthStats.singleMean = nanmean(spikeP75Width(singleSpikes));
spikeP75WidthStats.goodMean = nanmean(spikeP75Width(goodSpikes));
spikeP75WidthStats.singlePeakMean = nanmean(spikeP75Width(singleP25Spikes));
spikeP75WidthStats.allStd = nanstd(spikeP75Width);
spikeP75WidthStats.singleStd = nanstd(spikeP75Width(singleSpikes));
spikeP75WidthStats.goodStd = nanstd(spikeP75Width(goodSpikes));
spikeP75WidthStats.singlePeakStd = nanstd(spikeP75Width(singleP50Spikes));
% % start2peak amp - commenting out for now
% spikeSPAmpStats.allMean    = nanmean(spikeSPAmp);
% spikeSPAmpStats.singleMean = nanmean(spikeSPAmp(singleSpikes));
% spikeSPAmpStats.goodMean   = nanmean(spikeSPAmp(goodSpikes));
% spikeSPAmpStats.allStd     = nanstd(spikeSPAmp);
% spikeSPAmpStats.singleStd  = nanstd(spikeSPAmp(singleSpikes));
% spikeSPAmpStats.goodStd    = nanstd(spikeSPAmp(goodSpikes));
% % start2peak width
% spikeSPWidthStats.allMean    = nanmean(spikeSPWidth);
% spikeSPWidthStats.singleMean = nanmean(spikeSPWidth(singleSpikes));
% spikeSPWidthStats.goodMean   = nanmean(spikeSPWidth(goodSpikes));
% spikeSPWidthStats.allStd     = nanstd(spikeSPWidth);
% spikeSPWidthStats.singleStd  = nanstd(spikeSPWidth(singleSpikes));
% spikeSPWidthStats.goodStd    = nanstd(spikeSPWidth(goodSpikes));
% intracellular width
spikeIntraAmpStats.allMean    = nanmean(spikeIntraAmp);
spikeIntraAmpStats.singleMean = nanmean(spikeIntraAmp(singleSpikes));
spikeIntraAmpStats.goodMean   = nanmean(spikeIntraAmp(goodSpikes));
spikeIntraAmpStats.allStd     = nanstd(spikeIntraAmp);
spikeIntraAmpStats.singleStd  = nanstd(spikeIntraAmp(singleSpikes));
spikeIntraAmpStats.goodStd    = nanstd(spikeIntraAmp(goodSpikes));
% intracellular width
spikeIntraWidthStats.allMean    = nanmean(spikeIntraWidth);
spikeIntraWidthStats.singleMean = nanmean(spikeIntraWidth(singleSpikes));
spikeIntraWidthStats.goodMean   = nanmean(spikeIntraWidth(goodSpikes));
spikeIntraWidthStats.allStd     = nanstd(spikeIntraWidth);
spikeIntraWidthStats.singleStd  = nanstd(spikeIntraWidth(singleSpikes));
spikeIntraWidthStats.goodStd    = nanstd(spikeIntraWidth(goodSpikes));


%% MISCALLENOUS - firing rate, perfectPTWidth, perfectPTAMP, complex-spike index
% also characterize bursting...
% see http://jn.physiology.org/cgi/content/full/92/1/600
% they plotted spikeWidth vs spikePTWidth (spikePTWidth was the best
% variable for separation of spikes)
% mean firing rate... session variable is not available to this function,
% so just take the time between the first and last spikes. this is a close
% enough approximation to total duration - but remember that this is not
% precise and implement a more precise solution that uses the session
% variables to calculate exact duration of session
meanRate = length(spikeTimes) / (spikeTimes(end)-spikeTimes(1))

% perfectPTWidth: spikePTWidth appears to be the best measure of width,
% however since it looks for the first minima, it can be seriously thrown
% off by cases where there are multiple spikes. To
% correct for this and to assign a single perfectPTWidth value, find the
% spikes where spikePTWidth == spikeOverallMinPTWidth. For these values it
% is almost guaranteed that there was only one spike and that the spike was
% complete within the 1ms window. Hence this gives a very good average for
% perfectPTWidth
perfectSpikes = find(abs(spikeTPWidth - spikeOverallMaxTPWidth) < 1); % perfect is when there is less than 1 bin diff between firstpeak and overallpeak width. Before I had these set to equal. But because the firstpeak is found using the spikeSMOOTH and overallPeak is found using the non-smooth version, small differences are easily possible
perfectTPWidth = nanmean(spikeTPWidth(perfectSpikes));
perfectTPAmp   = nanmean(spikeTPAmp(perfectSpikes));


%Get correlation
%dTimes = diff(spikeTimes);
%newISIs = zeros(1, length(dTimes)+1);
%newISIs(2:end) = dTimes;
%goodISIs = find(newISIs(1:length(spikeWidth))/1e6 < 0.15);
%[isiWidthCorrAll pISIWidthAll] = corr(newISIs', spikeWidth');
%[isiWidthCorr150ms pISIWidth150ms] = corr(newISIs(goodISIs)', spikeWidth(goodISIs)');
%[widthAmpCorr pWidthAmpCorr] = corr(spikeWidth(singleSpikes)', spikeAmp(singleSpikes)');
%[isiAmpCorrAll pISIAmpCorrAll] = corr(newISIs', spikeAmp');
%[isiAmpCorr150ms pISIAmpCorr150ms] = corr(newISIs(goodISIs)', spikeAmp(goodISIs)');


%% CREATE WAVEFORMS
%allSingleWave = squeeze(ttData(:,modeMaxCh,spikeTTIndices(singleSpikes)));
%avgSingleWave = mean(ttData(:,modeMaxCh,spikeTTIndices(singleSpikes)),3);
%stdAvgSingleWave = std(allSingleWave,[],2);
allGoodWave = waveforms(:,goodSpikes);
avgGoodWave = mean(allGoodWave,2);
stdAvgGoodWave = std(allGoodWave,[],2);


%% SAVE VARIABLES
% 20090804 - moving save above figure generation so that data is saved in
% case something goes wrong in the figure generation
disp('saving....')
tic
disp('HERE')
%save('testSaveMat.mat','spikeT01Width')

	disp('here')
	save(cellPropSaveFile, ...
	    'spikeTimes', ...
	    'singleSpikes', 'dualSpikes', 'goodSpikes', 'perfectSpikes', 'nonNanSpikes', ...
	    'spikeMinIndex', 'spikeMaxIndex', 'spikeOverallMaxIndex', ... 'spikeStartMinIndex', ...
	    'spikeMinAmp', 'spikeMaxAmp', 'spikeT50Width','spikeT01Width','spikeT05Width', 'spikeT10Width','spikeT15Width','spikeT20Width','spikeT30Width','spikeT35Width','spikeT75Width','spikeTMinus05Width','spikeTMinus10Width','spikeTMinus15Width','spikeTMinus20Width', 'spikeTPAmp', 'spikeTPWidth', 'spikeOverallMaxTPAmp', 'spikeOverallMaxTPWidth', 'spikeIntraAmp', 'spikeIntraWidth','spikeIntraMinBin',...'spikeAvgTroughRiseDiff','spikeAvgPeakDecayDiff', ... 'spikeSPAmp', 'spikeSPWidth', 
	    'spikeP50Width', 'spikeP25Width', 'spikeP75Width', ...
	    'spikeMinAmpStats', 'spikeMaxAmpStats', 'spikeT50WidthStats', 'spikeTPAmpStats', 'spikeTPWidthStats', 'spikeOverallMaxTPAmpStats', 'spikeOverallMaxTPWidthStats', 'spikeIntraAmpStats', 'spikeIntraWidthStats', ... 'spikeSPAmpStats', 'spikeSPWidthStats',
	    'spikeP50WidthStats', 'spikeP25WidthStats', 'spikeP75WidthStats', ...
	    'spikeTP1090Width', 'spikeTP2080Width', ...
	    'spikeTP1090WidthStats', 'spikeTP2080WidthStats', ...
	    'perfectTPWidth', 'perfectTPAmp', ...
	    'avgGoodWave', 'stdAvgGoodWave', ...
	    'meanRate','waveforms');
%fds

if(0)
save(cellPropSaveFile, ...
    'spikeTimes', ...
    'singleSpikes', 'dualSpikes', 'goodSpikes', 'perfectSpikes', 'nonNanSpikes', ...
    'spikeMinIndex', 'spikeMaxIndex', 'spikeOverallMaxIndex', ... 'spikeStartMinIndex', ...
    'spikeMinAmp', 'spikeMaxAmp', 'spikeT50Width','spikeT01Width','spikeT05Width', 'spikeT10Width','spikeT15Width','spikeT20Width','spikeT30Width','spikeT35Width','spikeT75Width','spikeTMinus05Width','spikeTMinus10Width','spikeTMinus15Width','spikeTMinus20Width', 'spikeTPAmp', 'spikeTPWidth', 'spikeOverallMaxTPAmp', 'spikeOverallMaxTPWidth', 'spikeIntraAmp', 'spikeIntraWidth','spikeIntraMinBin',...'spikeAvgTroughRiseDiff','spikeAvgPeakDecayDiff', ... 'spikeSPAmp', 'spikeSPWidth', 
    'spikeTP1090Width', 'spikeTP2080Width', ...
    'meanRate','waveforms');
end
toc

if(plotFigures)
	%% FIGURES
	disp('plotting....')
	tic
	%% FIRST FIGURE
	%if cluster ~= 999 & cluster ~= 1000
	% one figure summing up the location of the peaks and the TP-widths
	binCenters = [0.25:0.5:48];
	fig = figure;
	fig_dim = [15 7.5];
	clf
	%set(fig, 'Color', [255 239 219]/255); % this sets the background color
	%set(fig, 'InvertHardCopy', 'off'); % this preserves the bg color when printing
	set(fig, 'PaperUnits', 'inches');
	set(fig, 'PaperSize', fig_dim);
	set(fig, 'PaperPosition', [0,0,(get(fig,'PaperSize'))])
	set(fig, 'visible', 'on')
	hold on; box on;
	% TPWidth Distribution
	subplot(1,3,1);
	hist(spikeOverallMaxTPWidth, binCenters);
	set(gca, 'fontsize', 15);
	xlabel('Trough-Peak Width (bins)');
	ylabel('Count');
	xlim([0 48]);
	% can put the following in the title if necessay: [escapedFolderName ',
	% % ' refTTCluster]
	meanString = sprintf('%.2f +- %.2f', spikeOverallMaxTPWidthStats.goodMean, spikeOverallMaxTPWidthStats.goodStd);
	titleString = sprintf('TPWidth = %s, full=%.2f, good=%.2f, perf=%.2f', meanString, length(nonNanSpikes)/length(spikeTimes), length(goodSpikes)/length(spikeTimes), length(perfectSpikes)/length(spikeTimes));
	title([titleString]);
	legendString = sprintf('PerfMean=%.2f', perfectTPWidth);
	legend(legendString);
	% TPAmp Distribution
	subplot(1,3,2);
	hist(spikeOverallMaxTPAmp, 50);
	set(gca, 'fontsize', 15);
	xlabel('Trough-Peak Amp (uV)');
	ylabel('Count');
	meanString = sprintf('%.2f +- %.2f', spikeOverallMaxTPAmpStats.goodMean, spikeOverallMaxTPAmpStats.goodStd);
	title(['PTAmp = ' meanString]);
	% TroughIndex Distribution
	subplot(1,3,3); hist(spikeMinIndex, binCenters);
	set(gca, 'fontsize', 15);
	xlabel('Trough Bin Number');
	ylabel('Count');
	xlim([0 numBins]);
	titleString = sprintf('PeakIndex, 6-16=%.2f', length(find(spikeMinIndex >= 6 & spikeMinIndex <= 16))/length(spikeTimes));
	title(titleString);
	% saveCellFigures(fig, 'spike_peak_trough', subjectName, folderName, TT, cluster, useChannel, channel);
	%cellPropPrintPrefix = sprintf('%s%g%s', saveDir, ch, currSuffix);
	%cellPropPrintPrefix = fullfile(saveDir,sprintf('%g%s', ch, currSuffix));
	cellPropPrintPrefix = fullfile(saveDir,sprintf('%g%s-%s', ch, currSuffix,subject));
	%cellPropPrintPrefix = fullfile(saveDir,sprintf('%g%s-%s-sz%d', ch, currSuffix,subject,szIdx));
	print('-djpeg100', '-r200', [cellPropPrintPrefix '_spike_trough_peak'])
	%print('-djpeg100', '-r200', [cellPropSaveFile '_spike_trough_peak'])
	%end

	%% SECOND FIGURE
	fig = figure;
	fig_dim = [12 8];
	clf
	%set(fig, 'Color', [255 239 219]/255); % this sets the background color
	%set(fig, 'InvertHardCopy', 'off'); % this preserves the bg color when printing
	set(fig, 'PaperUnits', 'inches');
	set(gcf, 'PaperSize', fig_dim);
	set(fig, 'PaperPosition', [0,0,(get(fig,'PaperSize'))])
	set(fig, 'visible', 'on')
	subplot(1,2,1);
	hold on; box on;
	set(gca, 'fontsize', 15);
	plot(avgGoodWave, 'b-', 'linewidth', 4);
	plot(-cumsum(avgGoodWave)/5, 'r', 'linewidth', 3);
	errorbar(avgGoodWave(:,1), stdAvgGoodWave);
	plot(avgGoodWave, 'b-', 'linewidth', 4);
	plot(-cumsum(avgGoodWave)/5, 'r', 'linewidth', 3);
	%    figureTitle = strcat(sprintf(['%s, %s: Max Channel Average Waveform. meanSpikeWidth: All=', num2str(meanSpikeWidth) ' Good=', num2str(meanSpikeWidthGood)], refTTCluster, folderName));
	%    escapedTitle = figureTitle;
	%    escapedTitle(regexp(escapedTitle, '_')) = ' ';
	%    title(escapedTitle);
	xlabel('Bin Number')
	ylabel('Amplitude (uV)')
	xlim([0 numBins]);
	legend('Extra', 'Intra/5');
	%print(gcf,'-dpng',[dataDir savePrefix 'Average_Waveform'])
	subplot(1,2,2);
	hold on; box on;
	set(gca, 'fontsize', 15);
	%plot(allGoodWave);
	plot(allGoodWave(:,1:13));
	plot(avgGoodWave, 'k-', 'linewidth', 10);
	plot(avgGoodWave, 'w-', 'linewidth', 2);
	%figureTitle = strcat(sprintf('%s, %s: Max Channel Peak 7 to 10 Waveforms and Average', refTTCluster, folderName));
	%escapedTitle = figureTitle;
	%escapedTitle(regexp(escapedTitle, '_')) = ' ';
	%title(escapedTitle);
	xlabel('Bin Number')
	ylabel('Amplitude (uV)')
	xlim([0 numBins]);
	% saveCellFigures(fig, 'spike_waveform', subjectName, folderName, TT, cluster, useChannel, channel);
	print('-djpeg100', '-r200', [cellPropPrintPrefix '_spike_waveform'])
	%print('-djpeg100', '-r200', [cellPropSaveFile '_spike_waveform'])
	toc
end

%% Load the data again so that it can be returned
%fds
clear spike* % free up some memory before loading the file
cellProp = load(cellPropSaveFile);

disp('Single run time:')
toc
return;
%% ADDITIONAL OLD ISI-WIDTH-AMP-CORRELATION FIGURES
%     figure; %Spikewidth versus Time
%     box on;
%     hold on;
%     set(gca, 'fontsize', 10);
%     plot((spikeTimes(dualSpikes)-spikeTimes(1))/1e6, spikeWidth(dualSpikes), 'r.')
%     plot((spikeTimes - spikeTimes(1))/1e6, spikeWidth, 'b.')
%     plot((spikeTimes(dualSpikes)-spikeTimes(1))/1e6, spikeWidth(dualSpikes), 'r.')
%     legend('Dual Spikes','location','best','orientation','vertical');
%     xlabel('Time (s)');
%     ylabel('Spike Width (bins)')
%     figureTitle = strcat(sprintf('%s, %s: Spike Width versus Time', refTTCluster, folderName));
%     escapedTitle = figureTitle;
%     escapedTitle(regexp(escapedTitle, '_')) = ' ';
%     title(escapedTitle);
%     print(gcf,'-dpng',[dataDir savePrefix 'SpikeWidth_Versus_Time'])
% 
%     figure; %SpikeAmp versus Time
%     box on;
%     hold on;
%     set(gca, 'fontsize', 10);
%     plot((spikeTimes(dualSpikes)-spikeTimes(1))/1e6, spikeAmp(dualSpikes)*1e6, 'r.')
%     plot((spikeTimes - spikeTimes(1))/1e6, spikeAmp*1e6, 'b.')
%     plot((spikeTimes(dualSpikes)-spikeTimes(1))/1e6, spikeAmp(dualSpikes)*1e6, 'r.')
%     legend('Dual Spikes','location','best','orientation','vertical');
%     xlabel('Time (s)');
%     ylabel('Spike Amplitude (microVolts)')
%     figureTitle = strcat(sprintf('%s, %s: Spike Amplitude versus Time', refTTCluster, folderName));
%     escapedTitle = figureTitle;
%     escapedTitle(regexp(escapedTitle, '_')) = ' ';
%     title(escapedTitle);
%     print(gcf,'-dpng',[dataDir savePrefix 'SpikeAmp_Versus_Time'])
% 
%     figure; %SpikeAmp versus SpikeWidth
%     box on;
%     hold on;
%     set(gca, 'fontsize', 10);
%     plot(spikeWidth(dualSpikes), spikeAmp(dualSpikes)*1e6, 'r.')
%     plot(spikeWidth, spikeAmp*1e6, 'b.')
%     plot(spikeWidth(dualSpikes), spikeAmp(dualSpikes)*1e6, 'r.')
%     legend('Dual Spikes','location','best','orientation','vertical');
%     xlabel('Spike Width (bins)');
%     ylabel('Spike Amplitude (microVolts)')
%     figureTitle = strcat(sprintf(['%s, %s: Spike Amplitude versus Spike Width. Corr = ' num2str(widthAmpCorr)], refTTCluster, folderName));
%     escapedTitle = figureTitle;
%     escapedTitle(regexp(escapedTitle, '_')) = ' ';
%     title(escapedTitle);
%     print(gcf,'-dpng',[dataDir savePrefix 'SpikeAmp_Versus_SpikeWidth'])
% 
%     figure; %Spikewidth versus ISI
%     box on;
%     hold on;
%     set(gca, 'fontsize', 10);
%     plot(newISIs(1:length(spikeWidth)), spikeWidth, 'b.')
%     xlabel('ISI (seconds)')
%     ylabel('Spike Width (bins)')
%     figureTitle = strcat(sprintf(['%s, %s: Spike Width versus ISI. rAll: ' num2str(isiWidthCorrAll) '. r150: ' num2str(isiWidthCorr150ms)], refTTCluster, folderName));
%     escapedTitle = figureTitle;
%     escapedTitle(regexp(escapedTitle, '_')) = ' ';
%     title(escapedTitle);
%     print(gcf,'-dpng',[dataDir savePrefix 'SpikeWidth_Versus_ISI'])
% 
%     figure; %SpikeAmp versus ISI
%     box on;
%     hold on;
%     set(gca, 'fontsize', 10);
%     plot(newISIs(1:length(spikeAmp)), spikeAmp*1e6, 'b.')
%     xlabel('ISI (seconds)')
%     ylabel('Spike Amplitude (microvolts)')
%     figureTitle = strcat(sprintf(['%s, %s: Spike Amplitude versus ISI. rAll: ' num2str(isiAmpCorrAll) '. r150: ' num2str(isiAmpCorr150ms)], refTTCluster, folderName));
%     escapedTitle = figureTitle;
%     escapedTitle(regexp(escapedTitle, '_')) = ' ';
%     title(escapedTitle);
%     print(gcf,'-dpng',[dataDir savePrefix 'SpikeAmp_Versus_ISI'])

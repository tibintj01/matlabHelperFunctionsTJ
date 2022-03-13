%% DETECT RIPPLES FROM HIPPOCAMPAL CSC
function rippleList = OmarrippleDetect(ripLFP,Fs)
% rippleCSC = 49;
% ripLFPFilename = strcat('../processed_data/', folderName, '/reshaped_csc', num2str(rippleCSC), '.mat');
% load(ripLFPFilename);
% eval(['ripLFP = csc' num2str(rippleCSC) ';']);
% Fs = CSC_SampleFrequency(1);

numBins = length(ripLFP); % The number of samples (i.e. points) in our dataset
lfpTimeAxis = (1:numBins) * (1/Fs); % This converts from bins to seconds

% now downsample
ds = 8;
ds = 1;
FsDS = Fs/ds; % downsample to 1 in 8
ripLFPDS = ripLFP(1:ds:end);
lfpTimeAxisDS = lfpTimeAxis(1:ds:end);
%cscFirstTS = CSC_Timestamp(1);

%% RIPPLE DETECTION PROPER
% ripple settings
peakSD = 6; % in non-smooth
shoulderSD = 1; % in smooth
smoothBins = 50; % 50 bins = 12.5 ms
lowFreq = 150;
highFreq = 250;
lfp = filterLFP(ripLFPDS, FsDS, lowFreq, highFreq, 2, 1);
smoothLFP = smooth(abs(lfp),smoothBins);
rippleIndices = find(abs(lfp) > peakSD & abs(lfp) < 80);

% figure;
% hold on;
% plot(lfpTimeAxisDS, zscore(ripLFPDS), 'g');
% plot(lfpTimeAxisDS, lfp)
% plot(lfpTimeAxisDS, smoothLFP, 'r');
% plot(lfpTimeAxisDS([1 end]), [shoulderSD shoulderSD], '--k');

%% GROUP CONTIGUOUS RIPPLES
% take the rippleIndices and run the shoulder detection on them first.
% for rippleIndices(1), calculate shoulders.
% now take the end shoulder timestamp as the cutoff
% keep incrementing rippleIndices until the timestamp of rippleIndices(i) > cutoff
% repeat process, this method is not dependent upon the minRippleSep variable below
numRipples = 0;
cutoff = 0;
rippleList = [];% zeros(1,1);
%figure(213)
for i=1:length(rippleIndices)
    if rippleIndices(i) < Fs % skip ripples that happen less than a second after the start of the session
        continue;
    end
    if rippleIndices(i) > cutoff
        numRipples = numRipples+1;
        % ripple detected, now go backwards in smoothLFP indices until value
        % is less than shoulderSD.
        startBin = rippleIndices(i);
        while smoothLFP(startBin) >= shoulderSD
            if startBin <= 1
                break;
            end
            startBin = startBin-1;
        end
        %startIndex = startIndex+1; % go forward to the next index where it is above the shoulderSD
        endBin = rippleIndices(i);
        while smoothLFP(endBin) >= shoulderSD
            if endBin >= length(smoothLFP)
                break;
            end
            endBin = endBin+1;
        end
        %endIndex = endIndex-1; % go back to the prev index where it is above the shoulderSD
        [peakValue peakBin] = max(abs(lfp(startBin:endBin))); % get the peak value and its location
        peakBin = peakBin + startBin - 1; % adjust the ripple index in terms of the main vector
        [peakSmoothValue peakSmoothBin] = max(smoothLFP(startBin:endBin)); % get the peak smooth value and its location
        peakSmoothBin = peakSmoothBin + startBin - 1; % adjust the ripple index in terms of the main vector
        duration = lfpTimeAxisDS(endBin) - lfpTimeAxisDS(startBin);
        
        
        rippleList(numRipples).startTime = lfpTimeAxisDS(startBin);
        rippleList(numRipples).peakTime  = lfpTimeAxisDS(peakBin);
        rippleList(numRipples).peakSmoothTime = lfpTimeAxisDS(peakSmoothBin);
        rippleList(numRipples).endTime = lfpTimeAxisDS(endBin);
        rippleList(numRipples).startBin = startBin;
        rippleList(numRipples).peakBin = peakBin;
        rippleList(numRipples).peakSmoothBin = peakSmoothBin;
        rippleList(numRipples).endBin = endBin;
        rippleList(numRipples).peakValue = peakValue;
        rippleList(numRipples).peakSmoothValue = peakSmoothValue;
        rippleList(numRipples).duration = duration;
        
        %rippleList(numRipples).rippleLFP=smoothLFP(startBin:endBin);
          %plot(lfp(startBin:endBin))
          %hold on
        
        cutoff = endBin;
    end
end

%% TEST PLOTS
%{
figure;
subplot(2,1,1)
histogram([rippleList.peakValue])
subplot(2,1,2)
histogram([rippleList.duration]);

figure

i = ceil(rand*length(rippleList))
% i = 53
timeAxis = (lfpTimeAxisDS(rippleList(i).peakBin-FsDS : rippleList(i).peakBin+FsDS) - lfpTimeAxisDS(rippleList(i).peakBin));
plot(timeAxis, abs(lfp(rippleList(i).peakBin-FsDS : rippleList(i).peakBin+FsDS)), 'r-', 'LineWidth', 2);
hold on;
plot(timeAxis, smoothLFP(rippleList(i).peakBin-FsDS:rippleList(i).peakBin+FsDS), 'c-', 'LineWidth', 4);
plot(timeAxis, zscore(ripLFPDS(rippleList(i).peakBin-FsDS:rippleList(i).peakBin+FsDS)), 'b-', 'LineWidth', 2);
plot(timeAxis, ones(1, length(timeAxis)), ':', 'LineWidth', 5, 'Color', [0 0 0]);;
startTimeOnAxis = (rippleList(i).startTime - lfpTimeAxisDS(rippleList(i).peakBin));
peakTimeOnAxis = (rippleList(i).peakTime - lfpTimeAxisDS(rippleList(i).peakBin));
endTimeOnAxis = (rippleList(i).endTime - lfpTimeAxisDS(rippleList(i).peakBin));
plot([startTimeOnAxis startTimeOnAxis], [-2 15], ':', 'LineWidth', 5, 'Color', [0 1 0]);
plot([peakTimeOnAxis peakTimeOnAxis], [-2 15], ':', 'LineWidth', 5, 'Color', [0.5 0.5 0.8]);
plot([endTimeOnAxis endTimeOnAxis], [-2 15], ':', 'LineWidth', 5, 'Color', [0 1 0]);
 xlim([-0.1 0.1])
 %}
end

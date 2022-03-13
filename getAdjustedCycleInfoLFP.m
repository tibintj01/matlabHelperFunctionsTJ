function [adjustedCycleInfo] = ...
    getAdjustedCycleInfoLFP(timeAxis,rawLFP, Fs, lowSlowFreq, highSlowFreq, lowFastFreq, highFastFreq,savePrefix,saveFilteredLFP)

PEAK_DETECT_THRESHOLD = 0.2;  % Added Dec 27, 2011 to make it more explicit that I am using a threshold of 0.2

if ~exist('saveFilteredLFP', 'var')
    saveFilteredLFP = 0;
end

if saveFilteredLFP == 1
    filteredLFPFilename = [savePrefix '_' num2str(lowSlowFreq) '-' num2str(highSlowFreq) '_adjustedBy_'  num2str(lowFastFreq) '-' num2str(highFastFreq) '_filtered_LFPs.mat'];
    if exist(filteredLFPFilename, 'file') ~= 2
        slowLfp = filterLFP(rawLFP, Fs, lowSlowFreq, highSlowFreq, 2, 1); % THIS IS A ZSCORED LFP
        slowLfp = round(slowLfp*1000)/1000; % round to 1/1000th of a Zscore. This seems to make almost no difference to the LFP. This saves File size space (3 fold reduction)
        
        fastLfp = filterLFP(rawLFP, Fs, lowFastFreq, highFastFreq, 2, 1); % THIS IS A ZSCORED LFP
        fastLfp = round(fastLfp*1000)/1000; % round to 1/1000th of a Zscore. This seems to make almost no difference to the LFP. This saves File size space (3 fold reduction)
        
        save(filteredLFPFilename, 'slowLfp','fastLfp');
        clear lfp
    end
end

filename = [savePrefix '_' num2str(lowSlowFreq) '-' num2str(highSlowFreq) '_minmax.mat'];

if (true || exist(filename, 'file') ~= 2)
    % filter and zscore the lfp
    if saveFilteredLFP == 1
        load(filteredLFPFilename);
    else
        slowLfp = filterLFP(rawLFP, Fs, lowSlowFreq, highSlowFreq, 2, 1);
        fastLfp = filterLFP(rawLFP, Fs, lowFastFreq, highFastFreq, 2, 1); % THIS IS A ZSCORED LFP
    end
    
    blockSize = Fs*60; % 1 min
    numBlocks = floor(length(slowLfp) / blockSize)
    numMaxPoints = 0;
    numMinPoints = 0;
    cycleBlockTransitions = 0;
    for i=1:numBlocks
        startIndex = (i-1)*blockSize+1;
        endIndex = i*blockSize;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %returns table of [idxes(:), vals(:)] for detected maxes and mins
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [currSlowLFPMaxPoints currSlowLFPMinPoints] = peakdet(slowLfp(startIndex:endIndex), PEAK_DETECT_THRESHOLD); 
        
        %dead time adjustment Tibin 11/14/2020
        if(isempty(currSlowLFPMaxPoints) || isempty(currSlowLFPMinPoints))
            continue
        end
        % adjust the indices so that the first index is a min and the last
        % is a max
        if currSlowLFPMinPoints(1,1) < currSlowLFPMaxPoints(1,1) % first column is index
            currSlowLFPMinPoints = currSlowLFPMinPoints(1:length(currSlowLFPMaxPoints),:); %make mins length match max length
        else
            currSlowLFPMaxPoints = currSlowLFPMaxPoints(2:end, :); %take out first max
            currSlowLFPMinPoints = currSlowLFPMinPoints(1:length(currSlowLFPMaxPoints),:); %make mins match max length
        end
        currLFPTimestamps = timeAxis(startIndex:endIndex);
        currNumMaxPoints = length(currSlowLFPMaxPoints);
        currNumMinPoints = length(currSlowLFPMinPoints);
        
        % Omar Addition - 20090624
        % There were 2 minor bugs here:
        % 1. The currRipMaxPoints being stored in RipMaxPoints were not the
        % real indices across the whole LFP. To fix this add the
        % startIndex-1 to each point
        cycleMaxPoints(numMaxPoints+1:numMaxPoints+currNumMaxPoints) = currSlowLFPMaxPoints(:,1) + startIndex - 1;
        cycleMinPoints(numMinPoints+1:numMinPoints+currNumMinPoints) = currSlowLFPMinPoints(:,1) + startIndex - 1;
        cycleMaxAmp(numMaxPoints+1:numMaxPoints+currNumMaxPoints) = currSlowLFPMaxPoints(:,2);
        cycleMinAmp(numMinPoints+1:numMinPoints+currNumMinPoints) = currSlowLFPMinPoints(:,2);
        cycleMaxTimes(numMaxPoints+1:numMaxPoints+currNumMaxPoints) = currLFPTimestamps(currSlowLFPMaxPoints(:,1));
        cycleMinTimes(numMinPoints+1:numMinPoints+currNumMinPoints) = currLFPTimestamps(currSlowLFPMinPoints(:,1));

        numMaxPoints = numMaxPoints + currNumMaxPoints;
        numMinPoints = numMinPoints + currNumMinPoints;
        cycleBlockTransitions(i) = numMaxPoints;
        %pause(0.001);
    end
    startIndex = numBlocks*blockSize + 1;
    endIndex = length(slowLfp);
    % only process if there is at least 1 seconds worth of data remaining
    if (endIndex - startIndex) > Fs
        [currSlowLFPMaxPoints currSlowLFPMinPoints] = peakdet(slowLfp(startIndex:endIndex), PEAK_DETECT_THRESHOLD);
        % adjust the indices so that the first index is a min and the last
        % is a max
        if currSlowLFPMinPoints(1,1) < currSlowLFPMaxPoints(1,1) % this is good
            currSlowLFPMinPoints = currSlowLFPMinPoints(1:length(currSlowLFPMaxPoints),:);
        else
            currSlowLFPMaxPoints = currSlowLFPMaxPoints(2:end, :);
            currSlowLFPMinPoints = currSlowLFPMinPoints(1:size(currSlowLFPMaxPoints,1),:);
        end
        currLFPTimestamps = timeAxis(startIndex:endIndex);
        currNumMaxPoints = length(currSlowLFPMaxPoints);
        currNumMinPoints = length(currSlowLFPMinPoints);
        cycleMaxPoints(numMaxPoints+1:numMaxPoints+currNumMaxPoints) = currSlowLFPMaxPoints(:,1) + startIndex - 1;
        cycleMinPoints(numMinPoints+1:numMinPoints+currNumMinPoints) = currSlowLFPMinPoints(:,1) + startIndex - 1;
        cycleMaxAmp(numMaxPoints+1:numMaxPoints+currNumMaxPoints) = currSlowLFPMaxPoints(:,2);
        cycleMinAmp(numMinPoints+1:numMinPoints+currNumMinPoints) = currSlowLFPMinPoints(:,2);
        cycleMaxTimes(numMaxPoints+1:numMaxPoints+currNumMaxPoints) = currLFPTimestamps(currSlowLFPMaxPoints(:,1));
        cycleMinTimes(numMinPoints+1:numMinPoints+currNumMinPoints) = currLFPTimestamps(currSlowLFPMinPoints(:,1));

        numMaxPoints = numMaxPoints + currNumMaxPoints;
        numMinPoints = numMinPoints + currNumMinPoints;
        cycleBlockTransitions(i) = numMaxPoints;
        %pause(0.001);
    end
    %% IMPORTANT
    % the above code guarantees that the first point is a minimum and the
    % last is a maximum

   
    %% DURATIONS
    cycleDurations = [diff(cycleMinTimes) 0]; % min-to-min, in SECONDS
    cycleMaxMaxDurations = [0 diff(cycleMaxTimes)];
    cycleMinMaxDurations = (cycleMaxTimes-cycleMinTimes);
    cycleMaxMinDurations = (cycleMinTimes(2:end)-cycleMaxTimes(1:end-1));
    
    %figure;
    %hold on;
    %plot(timeAxis, lfp, 'y');
    %plot(cycleMaxTimes, cycleMaxAmp, 'rx');
    %plot(cycleMinTimes, cycleMinAmp, 'gx');    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %ADJUST TO CLOSEST FAST RHYTHM PEAKS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    showPlots=0;
    
    peakLocalWindowSec=50e-3; 
    %peakLocalWindowSec=100e-3; 
     %peakLocalWindowSec=500e-3; 
    peakLocalHalfWindowInIdxes=round(peakLocalWindowSec*Fs/2); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %loop through each slow LFP peak and closest peak local to it in fastlfp
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     numMaxPts=length(cycleMaxPoints);
     lfpLength=length(slowLfp);
     
     disp('finding peaks in fast LFP........')
     tic
     [~,fastLfpPeakIdxes] = findpeaks(fastLfp);
     toc
     
     closestFastMaxIdxes=NaN(1,numMaxPts);

     
     disp('finding peaks in fast LFP corresponding to slow LFP peaks........')
     tic
    for mxi=1:numMaxPts
        
        currMaxIdx=cycleMaxPoints(mxi);
        
      [~,closestFastPkIdxIdx]=min(abs(currMaxIdx-fastLfpPeakIdxes));
      
      closestFastMaxIdxes(mxi)=fastLfpPeakIdxes(closestFastPkIdxIdx);
      
    
        
        
        %{
        if(showPlots)
                currMaxWindStartIdx=max(1,currMaxIdx-peakLocalHalfWindowInIdxes);
            currMaxWindEndIdx=min(lfpLength,currMaxIdx+peakLocalHalfWindowInIdxes);
        
         localFastLFP=fastLfp(currMaxWindStartIdx:currMaxWindEndIdx);
        %[~,peakLocalIdxes] = findpeaks(localFastLFP);
            figure
            localSlowLFP=slowLfp(currMaxWindStartIdx:currMaxWindEndIdx);
            plot(currMaxWindStartIdx:currMaxWindEndIdx,localFastLFP)
            hold on
            plot(currMaxWindStartIdx:currMaxWindEndIdx,localSlowLFP)
            
            vline(closestFastPeakIdx)
            %plot(peakLocalIdxes,peakValues,'k.','MarkerSize',30)
            close all
        end
      %}
        
    end
    
     toc
    
     closestFastMaxTimes=timeAxis(closestFastMaxIdxes);
    
    
    
    save(filename,'closestFastMaxTimes','closestFastMaxIdxes', 'cycleMaxAmp','cycleMinAmp','cycleDurations','cycleMaxMaxDurations','cycleMinMaxDurations','cycleMaxMinDurations',...
                                  'cycleMaxTimes','cycleMinTimes','cycleMaxPoints','cycleMinPoints', 'cycleBlockTransitions','lowSlowFreq', 'highSlowFreq', 'lowFastFreq', 'highFastFreq');
    adjustedCycleInfo = load(filename);
else
    adjustedCycleInfo = load(filename);
end

close all; clear all; clc
dataDir='/Users/tibinjohn/Downloads/hc-11/data';

sessionNames={'Achilles_10252013'	'Cicero_09012014' 'Gatsby_08022013' ...
'Achilles_11012013'	'Cicero_09102014'	'Gatsby_08282013'...
'Buddy_06272013' 'Cicero_09172014'};

saveDir='./hc11ThetaCycleData';

touchDir(saveDir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%loop through large binary files and save downsampled version during maze
%only
%Function to load various types of raw or processed data into matlab:
%{
bload.m
GetRawSpk.m
LoadBinary.m
LoadClu.m
LoadCluRes.m
LoadEvt.m
LoadFet.m
LoadSegs.m
LoadSpk.m
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fs=1250; %samples per second, see hc-11 data description 

lowFreq=6;
highFreq=12;
%minPeakRate=10; %Hz
minPeakRate=5; %Hz

totalNumCh=128;
loadChunkSize=16;
loadChunkSize=32;
%loadChunkSize=43;
%loadChunkSize=64;

numChunks=totalNumCh/loadChunkSize;

numSessions=length(sessionNames);



skipChannels=[];
for sessNum=1:numSessions
 
    
    currSessionName=sessionNames{sessNum};
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %see Channel Orders description pdf for skipped channels
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    skipChannels=[];
    if(contains(currSessionName,'Gatsby_08282013'))
        skipChannels=[35 45 47]-1; %xml starts counting channels at 0
    end
    if(contains(currSessionName,'Buddy_06272013'))
        skipChannels=[24 27 58]-1;%xml starts counting channels at 0
    end
    
    
    currDir=fullfile(dataDir,currSessionName);
    currDataFile=getFilePathsRegex(currDir,'*eeg');
    
    currDataInfoFile=getFilePathsRegex(currDir,'*mat');
    spikeInfo=load(currDataInfoFile);
    
    spikeStruct=spikeInfo.sessInfo.Spikes;
   positionStruct=spikeInfo.sessInfo.Position;
  
    pyrClusterIDs=spikeInfo.sessInfo.Spikes.PyrIDs;
    intClusterIDs=spikeInfo.sessInfo.Spikes.IntIDs;

    allClusterIDs=[pyrClusterIDs(:);intClusterIDs(:)];
    

    fileBase=fullfile(currDir,currSessionName);
    
    FsSpike=20000;
    
  
    currXMLdata=LoadXml(fileBase); 
    %T is time stamps in samples at FsSpike=20kHz
    %G is cluster number for each spike out of all clusters
    %Map gives original (spkGpNum,clusterNum) pair for each G(i) to get
    %waveforms
    %{
    disp('loading spike time and shape cluster data...')
    tic
    [TidxPerSpike,TotalClusterNumPerSpike,TotalCluNumToShankCluPair,~]=LoadCluRes(fileBase);
    toc
    %}

   
 
    
    mazeTaskTimeBounds=spikeInfo.sessInfo.Epochs.MazeEpoch;
    
    %inTaskSpikeTimeIdxes=spikeStruct.SpikeTimes>=mazeTaskTimeBounds(1) & spikeStruct.SpikeTimes<=mazeTaskTimeBounds(2);
   
    

     inTaskSpikeTimeIdxes=true(size(spikeStruct.SpikeTimes));
    allSpikeTimes=spikeStruct.SpikeTimes(inTaskSpikeTimeIdxes);
    allSpikeIDs=spikeStruct.SpikeIDs(inTaskSpikeTimeIdxes);
    
    unitIDnumbers=unique(spikeStruct.SpikeIDs);
    
 
    for gi=1:currXMLdata.nElecGps   
            shankIDStr=sprintf('%s_shank%d',currSessionName,gi);
            
            if(exist(sprintf('%s.txt',shankIDStr),'file'))
                disp('shank already processed, skipping...')
                continue
            end
            %clustersInGp=find(TotalCluNumToShankCluPair(:,2)==gi);
            
            %spikeIdxesInThisShank=TotalClusterNumPerSpike==totalClusterNum;
            
           %allSpikeTimesOnThisShank=TidxPerSpike(spikeIdxesInThisShank);
        
           %disp(sprintf('loading spike group %d waveform data',gi))
           tic
          numChannelsInGp=length(currXMLdata.SpkGrps(gi).Channels);
          numWaveformSamples=currXMLdata.SpkGrps(gi).nSamples;
        
          disp(sprintf('loading waveform shape for each suprathreshold event on shank %d...',gi))
         
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %missing data for skip channels
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          currChannels=[currXMLdata.SpkGrps(gi).Channels];
          currSkipChannels=intersect(skipChannels,currChannels);
          numChannelsInGp= numChannelsInGp-length(currSkipChannels); 
          
        [currGpWaveformData,numChannelsInGp]=LoadSpk([fileBase sprintf('.spk.%d',gi)],numChannelsInGp,numWaveformSamples);

        toc
        
        disp(sprintf('loading time and cluster ID of each suprathreshold event on shank %d...',gi))
        tic
        [TidxPerSpike,TotalClusterNumPerSpike,TotalCluNumToShankCluPair,~]=LoadCluRes(fileBase,[gi],[]);
        tPerSupraThresholdEvent=TidxPerSpike/FsSpike;
        
        numSpikesInShank=length(TidxPerSpike);
        totalCluNums= TotalCluNumToShankCluPair(:,1);
        realCluNums=TotalCluNumToShankCluPair(:,3);
        cluNumPerSpike=NaN(numSpikesInShank,1);
        for si=1:numSpikesInShank
            currSpikeTotCluNum=TotalClusterNumPerSpike(si);
            realCluNumIdx=find(totalCluNums==currSpikeTotCluNum);
            cluNumPerSpike(si)=realCluNums(realCluNumIdx);
        end
        cluNumPerSpike=cluNumPerSpike(:);
        
        actualClusterIDperSpike=100*gi+cluNumPerSpike; %100s place is shank num, rest in cluster# in shank
        toc
        disp('')
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %save LFP cycle info per channel into file
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %saveCycleInfoHC11

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %accumulate spike times for each cell and associated
        %channel for each into struct
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %spikeGroupNumPerSpikeTime=floor(allSpikeIDs/100);%see hc-11 data description
        %cluNumInGpPerSpikeTime=mod(allSpikeIDs,100);
     
        
        %allSpikeTimesInSpkGp=spikeStruct.SpikeTimes(spikeGroupNumPerSpikeTime==gi);
        %allSpikeIDsInSpkGp=spikeStruct.SpikeIDs(spikeGroupNumPerSpikeTime==gi);
        %thisGpThisCluSpikeTimes=TidxPerSpike(
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %which shank does current unit come from?
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        shankForCell=floor(allClusterIDs/100);%see hc-11 data description
        
        for i=1:length(unitIDnumbers)

            currIDnum=unitIDnumbers(i);
            saveFilePath=fullfile(saveDir,sprintf('%s_%d-%dTheta_Unit%dInfo.mat',currSessionName,lowFreq,highFreq,currIDnum));

            if(exist(saveFilePath,'file'))
                disp(sprintf('unit %d info file already exists, skipping...',currIDnum))
                continue
            end
            
            currShank=floor(currIDnum/100);

            %shankNum=floor(currIDnum/100); %see hc-11 data description
            %spikeGroupNum=spikeGroupNumPerSpikeTime(i);
            %spikeGroupNum=spikeGroupNumPerSpikeTime(i);
            if(currShank~=gi)
                continue %process only units in current spike group/shank
            end
            %clusterNumInShank=mod(currIDnum,100);%see hc-11 data description


            setCellTypeHC11

            thisUnitSpikeTimeIdx=(allSpikeIDs==currIDnum);
            unitSpikeTimes=allSpikeTimes(thisUnitSpikeTimeIdx);
            rateTimeBinWidth=0.05; %sec
            rateTimeBinWidth=1; %sec
            [rateTimeBinCenters,firingRateOverTime]=spikeTimesToFiringRate(unitSpikeTimes,rateTimeBinWidth);
            peakRate=max(firingRateOverTime)

            if(peakRate<minPeakRate)
                continue  %process only units with sufficiently high rates
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %find channel this unit comes from by which of 8 channels has
            %greatest amplitude
            %waveform data is nChannelsInShank x 32 x numSpikesInShank
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %collect idxes, waveforms corresponding to the spike times for this unit
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            thisUnitSpikeWaveformIdxes=(actualClusterIDperSpike==currIDnum);
            thisUnitSpikeWaveforms=currGpWaveformData(:,:,thisUnitSpikeWaveformIdxes);
            %thisUnitSpikeWaveforms=squeeze(thisUnitSpikeWaveforms);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %get average waveform shape for each channel in shank
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %figure; 
            avgWaveformAmpPerCh=zeros(numChannelsInGp,1);
            thisUnitSpikeAvgWaveforms=NaN(numChannelsInGp,size(currGpWaveformData,2));
            for ci=1:numChannelsInGp
                currChAvgWaveform=nanmean(thisUnitSpikeWaveforms(ci,:,:),3);
                thisUnitSpikeAvgWaveforms(ci,:)=currChAvgWaveform(:)';
                %plot(squeeze(thisUnitSpikeWaveforms(ci,:,1:10))+500*(ci-1))
                avgWaveformAmpPerCh(ci)=range(currChAvgWaveform);
                %plot(currChAvgWaveform+500*(ci-1))
                %hold on
            end
            %figure; plot(avgWaveformAmpPerCh)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %save channel number on this shank with greatest avg waveform amp
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [~,thisUnitCh]=max(avgWaveformAmpPerCh);
           
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %load LFP cycle data corresponding to this cell
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %currUnitCycleInfo=
            cycleInfoFileName=fullfile(saveDir,sprintf('%s_Ch%d_%d-%d_minmax.mat',currSessionName,thisUnitCh,lowFreq,highFreq));
            currUnitCycleInfo=load(cycleInfoFileName);
            
            minTimes=currUnitCycleInfo.cycleMinTimes;
                maxTimes=currUnitCycleInfo.cycleMaxTimes;
                
            
            numCycles=length(currUnitCycleInfo.cycleMaxTimes);
                [spikeAssignMax spikeAssignMin spikeAssignZero] = assignSpikesToCycles2017(unitSpikeTimes, minTimes, maxTimes);

                mode=2; %phase from max to max
                mode=1;%phase from startMin to endMin
                calcTimeOffsets=1;
                [spikePhase spikeOffsets spikeOffsetsEnd] = assignPhaseToSpikes(unitSpikeTimes, spikeAssignMax, spikeAssignMin, minTimes, maxTimes, mode, calcTimeOffsets, spikeAssignZero);
            
             
            unitSpikeThetaPhases=spikePhase(:);
            
            unitInfo.unitIDnum=currIDnum;
            unitInfo.unitCellType=cellType;
            unitInfo.unitShankNum=currShank;
            unitInfo.thisUnitSpikeAvgWaveforms=thisUnitSpikeAvgWaveforms;
            unitInfo.mazeTaskTimeBounds=mazeTaskTimeBounds;
            
            %unitInfo.spikeGroupNum=spikeGroupNum;
            unitInfo.thisUnitCh=thisUnitCh;
            unitInfo.cycleInfoFileName=cycleInfoFileName;
            unitInfo.positionPerTime=positionStruct;

            unitInfo.spikeTimes=unitSpikeTimes;
            unitInfo.spikePhases=unitSpikeThetaPhases;


            unitInfo.firingRateOverTime=firingRateOverTime;
            unitInfo.rateTimeBinCenters=rateTimeBinCenters;
            unitInfo.peakRate=peakRate;
            unitInfo.sessionName=currSessionName;
            
            save(saveFilePath,'unitInfo');
        end %unit ID loop

         shankFileID = fopen(sprintf('%s.txt',shankIDStr),'w');
        fprintf(shankFileID,sprintf('finished %s', shankIDStr))
        fclose(shankFileID)
    end%shank/spike group number
   
end%session loop

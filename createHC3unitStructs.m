close all; clear all; clc
dataDir='/Users/tibinjohn/hc-3';

metaDataFilePath='/Users/tibinjohn/hc-3/docs/hc3-metadata-tables/hc3-tables.db';
processedDataDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData';

channelPosMetaDataDir='/Users/tibinjohn/hc-3/docs/channelorder';

saveDir='./hc3ProcessedData/unitSpikeInfoStructs';
imageSaveDir='./hc3Images/unitLocationVsRippleLayer';
comboImgSaveDir='./hc3Images/placePhaseDepthComboImages';

saveLFPcycleInfoDir='./hc3ProcessedData/lfpCycleInfos';

savePositionInfoDir='./hc3ProcessedData/positionTimeInfos';

saveDirPlaceFieldAndPhaseImg='./hc3Images/placeFieldsAndPhaseImgs';

touchDir(saveDir)
touchDir(imageSaveDir)
touchDir(saveLFPcycleInfoDir)
touchDir(savePositionInfoDir)
touchDir(saveDirPlaceFieldAndPhaseImg)
touchDir(comboImgSaveDir)


NO_OVERWRITE=1;
%NO_OVERWRITE=0;
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
Fs=1250; %samples per second, see hc-3 data description 

lowFreq=6;
highFreq=12;
%minPeakRate=10; %Hz
minPeakRate=5; %Hz

tableColumnNames={'topdir','session','animal','duration'};

sqlQuery=['"SELECT distinct s.topdir, s.session, e.animal, s.duration ' ...
          'FROM cell c, session s, epos e, file f ' ...
          'WHERE c.topdir=s.topdir AND c.topdir=e.topdir AND s.session=f.session ' ...
          'AND c.region=''CA1'' '...
          'AND (s.behavior = ''linear'' OR s.behavior = ''linearOne'' OR s.behavior = ''linearTwo'') '...
          'ORDER BY s.duration desc;"'
]


%dbObj=sqlite3(metaDataFilePath,sqlQuery);
system(sprintf('sqlite3 %s %s > sqlQueries/sqlQueryResults.csv',metaDataFilePath,sqlQuery))

selectedSessions=readtable('sqlQueries/sqlQueryResults.csv')
selectedSessions.Properties.VariableNames=tableColumnNames;
%%
%numSessions=length(sessionNames);

numSessions=size(selectedSessions,1);

skipChannels=[];
for sessNum=1:numSessions
 
    currTopDir=selectedSessions{sessNum,1}{1};
    currSessionName=selectedSessions{sessNum,2}{1};
    
    if(NO_OVERWRITE && exist(sprintf('%s.txt',currSessionName),'file'))
        continue
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %get cell info from metadata database for this session
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    sessionTableColumnNames={'id','topdir', 'animal','ele','clu','region','fireRate','cellType'};
    sqlQueryForSessionPlaceholder=['"SELECT distinct c.id, c.topdir, c.animal, c.ele, c.clu, c.region, c.fireRate, c.cellType ' ...
          'FROM cell c, session s, epos e ' ...
          'WHERE c.topdir=''%s'' AND s.session=''%s'' AND c.cellType=''p'' ' ...
          'AND c.region=''CA1'' '...
          'AND (s.behavior = ''linear'' OR s.behavior = ''linearOne'' OR s.behavior = ''linearTwo'') '...
          'ORDER BY s.duration desc;"'
    ]

    sqlQueryForSession=sprintf(sqlQueryForSessionPlaceholder,currTopDir,currSessionName);
    
    system(sprintf('sqlite3 %s %s > sqlQueries/%s_sqlQueryResults.csv',metaDataFilePath,sqlQueryForSession,currSessionName))

    thisSessionTable=readtable(sprintf('sqlQueries/%s_sqlQueryResults.csv',currSessionName));
    thisSessionTable.Properties.VariableNames=sessionTableColumnNames;
    %thisSessionCellInfoTable=[];
    
    currDir=fullfile(dataDir,currTopDir,currTopDir,currSessionName);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %save position info for this session
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    posXYperTimeStepFileName=sprintf('%s.whl',currSessionName);
    posXYperTimeStepFilePath=fullfile(currDir,posXYperTimeStepFileName);
    if(~exist(posXYperTimeStepFilePath,'file'))
        continue
    end
    posXYperTimeStep=load(posXYperTimeStepFilePath);
    posXYperTimeStep(posXYperTimeStep==-1)=NaN;
    %maxDispIdx=15000;
    %{
    figure; 
    for maxDispIdx=2:150000
        if(~isnan(posXYperTimeStep(maxDispIdx,1)))
            plot(posXYperTimeStep((maxDispIdx-1):maxDispIdx,1), posXYperTimeStep((maxDispIdx-1):maxDispIdx,2),'b')
            hold on
            plot(posXYperTimeStep((maxDispIdx-1):maxDispIdx,3), posXYperTimeStep((maxDispIdx-1):maxDispIdx,4),'r')
            drawnow
        end
    end
    %}
    posXYperTimeStep=posXYperTimeStep(:,1:2); %use only head LED
    
    posFs=39.06; %samples per second
    posTimeStep=1/posFs; %seconds
    
    speedXPerTimeStep=movingslope(posXYperTimeStep(:,1),21,1); %~0.5 sec wind, cm per sample
    speedYPerTimeStep=movingslope(posXYperTimeStep(:,2),21,1); %~0.5 sec wind
    
    speedXPerTimeStep=speedXPerTimeStep*posFs; %cm per sec
    speedYPerTimeStep=speedYPerTimeStep*posFs; %cm per sec
    
    absSpeedPerTimeStep=sqrt(speedXPerTimeStep.^2 + speedYPerTimeStep.^2);
    
    %figure; histogram(absSpeedPerTimeStep)
    
    minSpeedDisp=40;%cm/sec
    
    highSpeedIdxes=absSpeedPerTimeStep>minSpeedDisp;
    
    [m,b,R]=getLinearFit(posXYperTimeStep(highSpeedIdxes,1),posXYperTimeStep(highSpeedIdxes,2))
    
    fitLineXmin=min(posXYperTimeStep(highSpeedIdxes,1));
    fitLineXmax=max(posXYperTimeStep(highSpeedIdxes,1));
    fitLineYmin=m*fitLineXmin+b;
    fitLineYmax=m*fitLineXmax+b;
    
    rightwardVectorAlongTrack=[(fitLineXmax-fitLineXmin); (fitLineYmax-fitLineYmin)];
    
    posAlongTrackPerTime=projectDataOntoVectorSpan(posXYperTimeStep',rightwardVectorAlongTrack);
    
    trackLengthCm=250; %cm
    trackLengthM=trackLengthCm/100; %cm
    maxTimeSec=length(posAlongTrackPerTime)/posFs; %sec
    
    posAlongTrackPerTimeCm=scaledata(posAlongTrackPerTime,0,trackLengthCm);
    
    posAlongTrackPerTimeM=posAlongTrackPerTimeCm/100;
    
    positionTimeAxisSec=(posTimeStep/2):posTimeStep:(maxTimeSec-posTimeStep/2);
    
    estSpeedWindSec=0.3; %sec
    estSpeedWindSec=1; %sec
    estSpeedWindIdx=round(estSpeedWindSec*posFs);
    %speedAlongTrackPerTime=movingslope(posAlongTrackPerTime,estSpeedWindIdx,1); %cm per time step
        speedAlongTrackPerTime=movingslope(posAlongTrackPerTimeCm,estSpeedWindIdx,2); %cm per time step

    speedAlongTrackPerTimeCmSec=speedAlongTrackPerTime*posFs;
    speedAlongTrackPerTimeMSec=speedAlongTrackPerTimeCmSec/100;
    
    
    %{
    figure; 
    yyaxis left
    plot(positionTimeAxis,posAlongTrackPerTime,'b')
    yyaxis right
    plot(positionTimeAxis,speedAlongTrackPerTimeCmSec,'r')
    %plot(posXYperTimeStep(highSpeedIdxes,1), posXYperTimeStep(highSpeedIdxes,2),'k.')
    hold on
    %plot([fitLineXmin fitLineXmax],[fitLineYmin fitLineYmax],'r-','LineWidth',4)
    title(currSessionName)
    daspect([1 1 1])
    %}
    positionInfoForSessionSavePath=fullfile(savePositionInfoDir,sprintf('%s_%s_PositionPerTimeInfo.mat',currTopDir, currSessionName));
    save(positionInfoForSessionSavePath,'positionTimeAxisSec', 'posAlongTrackPerTimeCm', 'posAlongTrackPerTimeM', 'speedAlongTrackPerTimeCmSec','speedAlongTrackPerTimeMSec', 'trackLengthCm','trackLengthM','estSpeedWindSec','posFs')
    %continue
    %save(savePositionInfoDir)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %FIND/READ PYRAMIDAL LAYER FOR TRUE THETA TIMING REFERENCE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [chPosMetaDataThisSessPath]=getChPosMetaDataFileForSess(channelPosMetaDataDir,currSessionName);
    
    maxRippleFileNameBase=getFileNameFromPath(chPosMetaDataThisSessPath);
    maxRippleFileNameBase=maxRippleFileNameBase(1:(end-4));
    pyrLayerChPerShank=load(chPosMetaDataThisSessPath);
    
    %{
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
    %}
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %LOAD LFPs for this session
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    %currDataFile=getFilePathsRegex(currDir,'*eeg');
    currEEGfilePath=getFilePathsRegex(currDir,'*eeg');
    fileBase=fullfile(currDir,currSessionName);
    
    currXMLdata=LoadXml(fileBase); 
    
    %FsSpike=20000;
    
    FsSpike=currXMLdata.SampleRate;
    
    totalNumCh=currXMLdata.nChannels;
    
    numElectrodesThisSession=currXMLdata.nElecGps;
    
    %totalNumCh=128;
    loadChunkSize=32;
    numChunks=ceil(totalNumCh/loadChunkSize);

    saveCycleInfoHC3
    %{
    saveRippleInfoHC3
    continue
    %}
    
    %{
    thisSessionRippleInfoPaths=getFilePathsRegex(processedDataDir,sprintf('%s_Ch*RippleInfo.mat',currSessionName));
    
    figure(87)
    rippleAmpEdges=linspace(0,30,50+1);
    rippleAmpBins=edgesToBins(rippleAmpEdges);
    meanRippleAmpPerCh=NaN(length(thisSessionRippleInfoPaths),1);
    for ci=1:length(thisSessionRippleInfoPaths)
        currRippleData=load(thisSessionRippleInfoPaths{ci});
        currRippleAmps=[currRippleData.rippleInfo.peakValue];
        meanRippleAmpPerCh(ci)=nanmedian(currRippleAmps);
        %{
        [N,c]=histcounts(currRippleAmps,rippleAmpEdges);
        plot(rippleAmpBins,N/sum(N))
        hold on
        %}
    end
    plot(meanRippleAmpPerCh,'ko')
    disp('')
    %}

    
 
    %}
    
    %mazeTaskTimeBounds=spikeInfo.sessInfo.Epochs.MazeEpoch;
    
    %inTaskSpikeTimeIdxes=spikeStruct.SpikeTimes>=mazeTaskTimeBounds(1) & spikeStruct.SpikeTimes<=mazeTaskTimeBounds(2);
   
    
%{
     inTaskSpikeTimeIdxes=true(size(spikeStruct.SpikeTimes));
    allSpikeTimes=spikeStruct.SpikeTimes(inTaskSpikeTimeIdxes);
    allSpikeIDs=spikeStruct.SpikeIDs(inTaskSpikeTimeIdxes);
    
    unitIDnumbers=unique(spikeStruct.SpikeIDs);
    %}
 
       %T is time stamps in samples at FsSpike=20kHz
    %G is cluster number for each spike out of all clusters
    %Map gives original (spkGpNum,clusterNum) pair for each G(i) to get
    %waveforms
    
    %{
    disp('loading all spike times, waveforms, and shape cluster data...')
    tic
        [TidxPerSpike,TotalClusterNumPerSpike,TotalCluNumToShankCluPair,~]=LoadCluRes(fileBase);
    %}

        shanks={currXMLdata.SpkGrps.Channels};
        nSamplesPerShank={currXMLdata.SpkGrps.nSamples};

        %waveformDataFile=sprintf('clusterAndWaveformDataSession%s.mat',currSessionName);
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %load all waveform info for this session
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tic
        for ei=1:numElectrodesThisSession
             disp(sprintf('loading time and cluster ID of each suprathreshold event on shank %d...',ei))
            
            [TidxPerSpike,TotalClusterNumPerSpike,TotalCluNumToShankCluPair,~]=LoadCluRes(fileBase,[ei],[]);
            
            if(isempty(TidxPerSpike))
                continue
            end
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

            actualClusterIDperSpike=100*ei+cluNumPerSpike; %100s place is shank num, rest in cluster# in shank
            %toc
        
            numChannelsInGp=length(shanks{ei});
            numWaveformSamples=nSamplesPerShank{ei};

            [currGpWaveformData,numChannelsInGp]=LoadSpk([fileBase sprintf('.spk.%d',ei)],numChannelsInGp,numWaveformSamples);

            allWaveformsPerShank{ei}=currGpWaveformData;
            allShankCluNumIDsPerShank{ei}=actualClusterIDperSpike;
            TidxPerSpikePerShank{ei}=TidxPerSpike;
            numChannelsInGpPerShank{ei}=numChannelsInGp;
        end
        
        
    
    toc
    disp('')
    %%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %loop through selected pyr CA1 cells to get spike times, closest ch
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for ci=1:size(thisSessionTable,1)
        currCellShankNum=thisSessionTable{ci,4};
        currCellCluNumInShank=thisSessionTable{ci,5};
        
        currCellShankCluID=currCellShankNum*100+currCellCluNumInShank;
        
        saveFilePath=fullfile(saveDir,sprintf('%s_Unit%dInfo.mat',currSessionName,currCellShankCluID));
    
        if(NO_OVERWRITE && exist(saveFilePath,'file'))
            continue
        end
        
        currShankCluNumIDs=allShankCluNumIDsPerShank{currCellShankNum};
        currGpWaveformData=allWaveformsPerShank{currCellShankNum};
        TidxPerSpike=TidxPerSpikePerShank{currCellShankNum};
        
        getPyrLayerIdx=find(currCellShankNum==pyrLayerChPerShank(:,1));
        maxRippleChThisShank=pyrLayerChPerShank(getPyrLayerIdx,2);
        
        
        chsThisShank=[shanks{currCellShankNum}]+1;
        
        if(length(chsThisShank)<8)
            disp('')
        end
        
        
       numChPerShank=max([numChannelsInGpPerShank{:}]);
       % maxRippleLevelThisShank=mod(maxRippleChThisShank,numChannelsInGpPerShank{currCellShankNum});
       %maxRippleLevelThisShank=mod(maxRippleChThisShank,numChPerShank)+1;
       maxRippleLevelThisShank=find(chsThisShank==maxRippleChThisShank);
       
       %{
       if(maxRippleLevelThisShank==0)
           maxRippleLevelThisShank=numChPerShank;
       end
       %}

        
        %currCellTotalCluNum=find(currShankCluNumIDs==currCellShankCluID & TotalCluNumToShankCluPair(:,3)==currCellCluNumInShank);
    
        currCellSpikeTimeInIdxes=TidxPerSpike(currShankCluNumIDs==currCellShankCluID);
        currCellSpikeTimesSec=currCellSpikeTimeInIdxes/FsSpike;
        unitSpikeTimes=currCellSpikeTimesSec;
        
          posInfo=load(positionInfoForSessionSavePath);
        [positionBins,firingRatePerPositionRight, firingRatePerPositionLeft,spikeTimeDirectionAssignment] ...
               = getSpikeStatsPerPosition(unitSpikeTimes,posInfo.positionTimeAxisSec,posInfo.posAlongTrackPerTimeM,posInfo.speedAlongTrackPerTimeMSec)
        
        if(max([firingRatePerPositionRight(:);firingRatePerPositionLeft(:)])<1)
            disp('LESS THAN 1HZ spatial firing, skipping cell')
            continue
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %find channel this unit comes from by which of 8 channels has
        %greatest amplitude
        %waveform data is nChannelsInShank x 32 x numSpikesInShank
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %collect idxes, waveforms corresponding to the spike times for this unit
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            thisUnitSpikeWaveformIdxes=(currShankCluNumIDs==currCellShankCluID);
            thisUnitSpikeWaveforms=currGpWaveformData(:,:,thisUnitSpikeWaveformIdxes);
            %thisUnitSpikeWaveforms=squeeze(thisUnitSpikeWaveforms);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %get average waveform shape for each channel in shank
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            figure; 
            numChannelsInGp=numChannelsInGpPerShank{currCellShankNum};
            
            avgWaveformAmpPerCh=zeros(numChannelsInGp,1);
            thisUnitSpikeAvgWaveforms=NaN(numChannelsInGp,size(currGpWaveformData,2));
            for ci=1:numChannelsInGp
                currChAvgWaveform=nanmean(thisUnitSpikeWaveforms(ci,:,:),3);
                thisUnitSpikeAvgWaveforms(ci,:)=currChAvgWaveform(:)';
                maxNumSpikesDisp=min(100,size(thisUnitSpikeWaveforms(ci,:,:),3));
                
                correctedCi=CorrectChaOrderCA1(maxRippleFileNameBase,currCellShankNum,ci)

                plot((squeeze(thisUnitSpikeWaveforms(ci,:,1:maxNumSpikesDisp))+500*(correctedCi))/500)
                avgWaveformAmpPerCh(ci)=range(currChAvgWaveform);
                
                %plot((squeeze(thisUnitSpikeWaveforms(ci,:,1:maxNumSpikesDisp))+500*(ci))/500)
           
                %plot((currChAvgWaveform+500*(ci-1))/500)
                hold on
            end
           xlim([1 size(thisUnitSpikeWaveforms,2)])
           ylim([0 9])
            
            xlabel('waveform sample num')
            ylabel('channel depth within shank')
            
            %figure; plot(avgWaveformAmpPerCh)

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %save channel number on this shank with greatest avg waveform amp
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [~,thisUnitChLevel]=max(avgWaveformAmpPerCh);
            thisUnitCOMch=getCenterOfMass(avgWaveformAmpPerCh);
            
            
            maxRippleLevelThisShank = CorrectChaOrderCA1(maxRippleFileNameBase,currCellShankNum,maxRippleLevelThisShank)
            thisUnitChLevel = CorrectChaOrderCA1(maxRippleFileNameBase,currCellShankNum,thisUnitChLevel)

            chLevelRelativeToPyrLayer=thisUnitChLevel-maxRippleLevelThisShank;

                        
            %midX=range(xlim)/2;
            thisUnitH=plot(xlim,[thisUnitChLevel thisUnitChLevel],'r-','LineWidth',5);
            %thisUnitH=plot(xlim,[thisUnitCOMch thisUnitCOMch],'r-','LineWidth',5);
            
            ripMaxLevelH=plot(xlim,[maxRippleLevelThisShank maxRippleLevelThisShank]+0.05,'k-','LineWidth',5);

            legend([thisUnitH ripMaxLevelH ],'this unit level (max spike amp)','pyramidal cell layer ch. (max ripple amp)','Location','northeast')
        
            title(sprintf('%s Cell %d depth profile across shank',currSessionName, currCellShankCluID))
            setFigFontTo(14)
            
            imageFileName=sprintf('%s_Cell%dDepthProfile.tif',currSessionName, currCellShankCluID);
            depthProfileImgPath=fullfile(imageSaveDir,imageFileName);
            saveas(gcf,depthProfileImgPath)
            close all
        
       
            cycleInfoFileName=fullfile(saveLFPcycleInfoDir,sprintf('%s_Ch%d_%d-%d_minmax.mat',currSessionName,maxRippleChThisShank,lowFreq,highFreq));
            currUnitCycleInfo=load(cycleInfoFileName);
            
            minTimes=currUnitCycleInfo.cycleMinTimes;
            maxTimes=currUnitCycleInfo.cycleMaxTimes;
                
            
            numCycles=length(currUnitCycleInfo.cycleMaxTimes);
            [spikeAssignMax spikeAssignMin spikeAssignZero] = assignSpikesToCycles2017(unitSpikeTimes, minTimes, maxTimes);

            mode=2; %phase from max to max
            %mode=1;%phase from startMin to endMin
            calcTimeOffsets=1;
            [spikePhase spikeOffsets spikeOffsetsEnd] = assignPhaseToSpikes(unitSpikeTimes, spikeAssignMax, spikeAssignMin, minTimes, maxTimes, mode, calcTimeOffsets, spikeAssignZero);

            unitSpikeThetaPhases=spikePhase(:);
            
            unitInfo.unitIDnum=currCellShankCluID;
            unitInfo.unitCellType='pyr';
            unitInfo.unitShankNum=currCellShankNum;
            %unitInfo.thisUnitSpikeAvgWaveforms=thisUnitSpikeAvgWaveforms;
            %unitInfo.mazeTaskTimeBounds=mazeTaskTimeBounds;
            
            %unitInfo.spikeGroupNum=spikeGroupNum;
            unitInfo.thisUnitChLevel=thisUnitChLevel;
            unitInfo.maxRippleLevelThisShank=maxRippleLevelThisShank;
            unitInfo.maxRippleChThisShank=maxRippleChThisShank;
            %unitInfo.chLevelRelativeToPyrLayer=chLevelRelativeToPyrLayer;
            %unitInfo.cycleInfoFileName=cycleInfoFileName;
            %unitInfo.positionPerTime=positionStruct;

            unitInfo.spikeTimes=unitSpikeTimes;
            unitInfo.spikePhases=unitSpikeThetaPhases;


            %unitInfo.firingRateOverTime=firingRateOverTime;
            %unitInfo.rateTimeBinCenters=rateTimeBinCenters;
            %unitInfo.peakRate=peakRate;
            unitInfo.sessionName=currSessionName;
            unitInfo.positionInfoForSessionSavePath=positionInfoForSessionSavePath;
            
            unitInfo.depthProfileImgPath=depthProfileImgPath;
            
          
           
           unitInfo.positionBins=positionBins;
           unitInfo.firingRatePerPositionRight=firingRatePerPositionRight;
           unitInfo.firingRatePerPositionLeft=firingRatePerPositionLeft;
           unitInfo.spikeTimeDirectionAssignment=spikeTimeDirectionAssignment;
           
           spikePositions=interp1(posInfo.positionTimeAxisSec,posInfo.posAlongTrackPerTimeM,unitSpikeTimes);
            
           rightSpikePositions=spikePositions(spikeTimeDirectionAssignment==1);
           rightSpikePhases=unitSpikeThetaPhases(spikeTimeDirectionAssignment==1);
           
           leftSpikePositions=spikePositions(spikeTimeDirectionAssignment==2);
           leftSpikePhases=unitSpikeThetaPhases(spikeTimeDirectionAssignment==2);
           
           unitInfo.spikePositions=spikePositions;
           
           unitInfo.leftSpikePositions=leftSpikePositions;
           unitInfo.leftSpikePhases=leftSpikePhases;
           
           unitInfo.rightSpikePositions=rightSpikePositions;
           unitInfo.rightSpikePhases=rightSpikePhases;
           
           unitInfo.speedPerTimeStep=speedAlongTrackPerTimeMSec;
           unitInfo.positionPerTimeStep=posAlongTrackPerTimeM;
           
           saveImgPathRightward=plotHC3PlaceFieldsAndThetaPhases(unitInfo,'betterRightward',saveDirPlaceFieldAndPhaseImg)
           saveImgPathLeftward=plotHC3PlaceFieldsAndThetaPhases(unitInfo,'betterLeftward',saveDirPlaceFieldAndPhaseImg)
           
           unitInfo.saveImgPathRightward=saveImgPathRightward;
           unitInfo.saveImgPathLeftward=saveImgPathLeftward;
           
           figure; 
           setTightSubplots
           subplot(2,2,1)
           try
            imshow(saveImgPathRightward)
           end
           title('rightward laps')
           subplot(2,2,3)
           try
            imshow(saveImgPathLeftward)
           end
           title('leftward laps')
           subplot(2,2,[2 4])
           try
            imshow(depthProfileImgPath)
           end
           uberTitle(sprintf('%s Cell %d (%s)',currSessionName,currCellShankCluID, unitInfo.unitCellType))
           setFigFontTo(14)
           maxFig
           
           comboImageFileName=sprintf('%s_Cell%dCombo.tif',currSessionName, currCellShankCluID);
           comboImageFilePath=fullfile(comboImgSaveDir,comboImageFileName);
           saveas(gcf,comboImageFilePath)
           close all
           
           unitInfo.comboImageFilePath=comboImageFilePath;
           
            save(saveFilePath,'unitInfo');
        
    end
    
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
        %shankForCell=floor(allClusterIDs/100);%see hc-11 data description
        
       
            
            
      

         sessionFileID = fopen(sprintf('%s.txt',currSessionName),'w');
        fprintf(sessionFileID,sprintf('finished %s', currSessionName))
        fclose(sessionFileID)
    %end%shank/spike group number
   
end%session loop



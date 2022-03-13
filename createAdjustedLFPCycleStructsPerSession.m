close all; clear all; clc
dataDir='/Users/tibinjohn/hc-3';

metaDataFilePath='/Users/tibinjohn/hc-3/docs/hc3-metadata-tables/hc3-tables.db';
processedDataDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData';

channelPosMetaDataDir='/Users/tibinjohn/hc-3/docs/channelorder';

saveDir='./hc3ProcessedData/unitSpikeInfoStructs';
imageSaveDir='./hc3Images/unitLocationVsRippleLayer';
comboImgSaveDir='./hc3Images/placePhaseDepthComboImages';

saveLFPadjustedCycleInfoDir='./hc3ProcessedData/lfpOmarAsymAdjustedCycleInfos';

saveRawLFPInfoDir='./hc3ProcessedData/rawLFPstructs';

saveDirPlaceFieldAndPhaseImg='./hc3Images/placeFieldsAndPhaseImgs';

touchDir(saveDir)
touchDir(imageSaveDir)
touchDir(saveLFPadjustedCycleInfoDir)
touchDir(saveRawLFPInfoDir)
touchDir(saveDirPlaceFieldAndPhaseImg)
touchDir(comboImgSaveDir)

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
    
    pyrLayerChannelNums=[];
    currTopDir=selectedSessions{sessNum,1}{1};
    currSessionName=selectedSessions{sessNum,2}{1};
    
   
    %lfpStructForSessionSavePath=fullfile(saveRawLFPInfoDir,sprintf('%s_%s_adjCycleLFPstruct.mat',currTopDir, currSessionName));
        lfpStructForSessionSavePath=fullfile(saveLFPadjustedCycleInfoDir,sprintf('%s_%s_adjCycleLFPstruct.mat',currTopDir, currSessionName));


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %FIND/READ PYRAMIDAL LAYER FOR TRUE THETA TIMING REFERENCE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [chPosMetaDataThisSessPath]=getChPosMetaDataFileForSess(channelPosMetaDataDir,currSessionName);
    if(isnan(chPosMetaDataThisSessPath))
        continue
    end
    
    maxRippleFileNameBase=getFileNameFromPath(chPosMetaDataThisSessPath);
    maxRippleFileNameBase=maxRippleFileNameBase(1:(end-4));
    pyrLayerChPerShank=load(chPosMetaDataThisSessPath);
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %LOAD LFPs for this session
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    currDir=fullfile(dataDir,currTopDir,currTopDir,currSessionName);
    
    currEEGfilePath=getFilePathsRegex(currDir,'*eeg');
    fileBase=fullfile(currDir,currSessionName);
    
    currXMLdata=LoadXml(fileBase); 

    FsSpike=currXMLdata.SampleRate;
    
    totalNumCh=currXMLdata.nChannels;
    
    numElectrodesThisSession=currXMLdata.nElecGps;
    
    %totalNumCh=128;
    loadChunkSize=32;
    numChunks=ceil(totalNumCh/loadChunkSize);
    
    pyrLayerChannelNums=pyrLayerChPerShank(:,2);

    %saveFullLFPinfoHC3
    %saveGammaCycleInfoHC3
    %saveAllChGammaCycleInfoHC3
    
    saveAdjustedCycleInfoHC3
         
   
end%session loop



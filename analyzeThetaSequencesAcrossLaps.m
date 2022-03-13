close all;  clc
%dataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles';

clearvars -except spikeDataPerField

tic
disp('loading spike data per field')
if(~exist('spikeDataPerField','var'))
    spikeDataPerField=load('allUnitInFieldGoodSpeedSpikePhasesAndTimes.mat');
end
toc

maxRunSpeed=1;


processedDataDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData';
unitDataDir=fullfile(processedDataDir,'unitSpikeInfoStructs');
dataDir=unitDataDir;

saveThetaSeqImgDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3Images/individualThetaSequencesSimilarSizeMinTrajLengthHalfWidth';
touchDir(saveThetaSeqImgDir)

posInfoDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData/positionTimeInfos';
perCycleDataDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData/perFieldPerCycleData_Dec16_2021';

thetaSequenceImgDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData/thetaSequenceImgDir';
touchDir(thetaSequenceImgDir)


%cycleDataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/';
showPlots=1;
%showPlots=0;



sqlDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/sqlQueries';

sqlSessionFiles=getRegexFileNames(sqlDir,'ec*csv');

numSessions=length(sqlSessionFiles);
sessionNames=cell(numSessions,1);

for si=1:numSessions
    sessionNames{si}=getSubstrBtwnSubstrs(sqlSessionFiles{si},'','_sql');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%loop through each session and load per cycle spike time data across units
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

startSi=1;
%startSi=76;

minNumFieldsInSeq=3;
minNumFieldsInSeq=4;
%minNumFieldsInSeq=2;
%startSi=16;
    
onlyCorrectOrder=1;
onlyCorrectOrder=0;

onlySimilarSizeFields=1;
%minAllowableFieldRatio=0.8;
minAllowableFieldRatio=0.7;

enforceMinTrajectoryLength=1;
%minTrajLengthWidthRatio=1/3;
minTrajLengthWidthRatio=1/2;

minNumOccurences=5;

setTightSubplots_Medium

for si=startSi:length(sessionNames) 
    currSessName=sessionNames{si};
    currSeqDataFilePath=fullfile(perCycleDataDir,sprintf('perFieldPerCycleTimingData%s.mat',currSessName));
    
    %missing seq files?
    if(exist(currSeqDataFilePath,'file'))
        currThetaSeqData=load(currSeqDataFilePath)
    else
        continue
    end
    
    if(isempty(currThetaSeqData.thetaSequenceInfoPerSequence))
        continue
    end
    
    thetaSequenceInfoPerSequence=currThetaSeqData.thetaSequenceInfoPerSequence;
    
    seqStrs=fieldnames(thetaSequenceInfoPerSequence);
    
    for ti=1:length(seqStrs)
        currSeqData=thetaSequenceInfoPerSequence.(seqStrs{ti});
        numOccurences=currSeqData.numOccurences;
        
        %if(numOccurences>=minNumOccurences)
            cycleNumPerOcc=currSeqData.cycleNumPerOcc;
            currThetaSeqSummaryImgPaths=[];
            for oi=1:length(cycleNumPerOcc)
                currCycleNum=cycleNumPerOcc{oi};
                seqImgPath=getFilePathsRegex(saveThetaSeqImgDir,sprintf('%s*Cycle%d_*',currSessName,currCycleNum));
                
                currThetaSeqSummaryImgPaths{oi}=seqImgPath;
               
            end
            
            %update theta seq info with summary image path array
            thetaSequenceInfoPerSequence.(seqStrs{ti}).currThetaSeqSummaryImgPaths=currThetaSeqSummaryImgPaths;
            
            %save(currSeqDataFilePath
            
            disp('')
            
        %end
    end
    
    %save updated theta seq info for all sequences, with summary image path array
    save(currSeqDataFilePath,'thetaSequenceInfoPerSequence','-append')
    
end
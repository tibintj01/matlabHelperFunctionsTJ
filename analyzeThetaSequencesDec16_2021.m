close all;  clc
%dataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles';

clearvars -except spikeDataPerField
%clear all

tic
disp('loading spike data per field')
if(~exist('spikeDataPerField','var'))
    spikeDataPerField=load('allUnitInFieldGoodSpeedSpikePhasesAndTimes.mat');
end
toc

maxRunSpeed=1;
maxRunSpeed=0.8;


processedDataDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData';
unitDataDir=fullfile(processedDataDir,'unitSpikeInfoStructs');
dataDir=unitDataDir;

%saveThetaSeqImgDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3Images/individualThetaSequencesAnyTrajLengthMinNumFields_4_WithFiringRate';
saveThetaSeqImgDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3Images/individualThetaSequencesAnyTrajLengthMinNumFields2';

touchDir(saveThetaSeqImgDir)

useOmarAsymAdjustment=1;

if(useOmarAsymAdjustment)
     adjustedCycleDataDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData/lfpOmarAsymAdjustedCycleInfos';
else
    adjustedCycleDataDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData/lfpAdjustedCycleInfos';
end

posInfoDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData/positionTimeInfos';
%perCycleDataDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData/perCycleUnitSeqData';
%perCycleDataDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData/perFieldPerCycleData';
perCycleDataDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData/perFieldPerCycleData_Dec16_2021';

thetaSequenceImgDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData/thetaSequenceImgDir';
touchDir(thetaSequenceImgDir)


%cycleDataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/';
showPlots=1;
showPlots=0;

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



minNumFieldsInSeq=2;
%minNumFieldsInSeq=3;
%minNumFieldsInSeq=4;
%minNumFieldsInSeq=5;
%minNumFieldsInSeq=2;
%startSi=16;
    
onlyCorrectOrder=1;
onlyCorrectOrder=0;

onlySimilarSizeFields=1;
minAllowableFieldRatio=0.8;
%minAllowableFieldRatio=0.7;
%minAllowableFieldRatio=0.75;
%minAllowableFieldRatio=0.9;
minAllowableFieldRatio=0.75;


%minAllowableFieldRatio=1/3;

enforceMinTrajectoryLength=1;
%minTrajLengthWidthRatio=1/3;
minTrajLengthWidthRatio=1/2;
minTrajLengthWidthRatio=1/3;
minTrajLengthWidthRatio=1/5;
minTrajLengthWidthRatio=1/8; %at least ~1 gamma cycle of center difference
minTrajLengthWidthRatio=0;

setTightSubplots_Medium

startSi=1;
%startSi=76;

for si=startSi:length(sessionNames) 
    
    thetaSequenceInfoPerSequence=[];
    
    si
    currSessName=sessionNames{si};
    
    %{
    if(~strcmp(currSessName,'ec013.626'))
        continue
    end
    %}
    currSeqDataFilePath=fullfile(perCycleDataDir,sprintf('perFieldPerCycleTimingData%s.mat',currSessName));
    
    currPosInfoFilePath=getFilePathsRegex(posInfoDir,sprintf('*%s*Position*.mat',currSessName));
    
    currSessUnitFilePaths=getFilePathsRegex(unitDataDir,sprintf('*%s*.mat',currSessName));

    if(isempty(currSessUnitFilePaths))
        continue
    end
    
    loadRepresentativeUnitFilesForSess
    
    
    if(~isempty(representativeLeftwardUnitInfo))
        leftwardLapStartTimes=representativeLeftwardUnitInfo.lapStartTimesPerDir.leftward;
        leftwardLapStopTimes=representativeLeftwardUnitInfo.lapStopTimesPerDir.leftward;
    end
    
    if(~isempty(representativeRightwardUnitInfo))
        rightwardLapStartTimes=representativeRightwardUnitInfo.lapStartTimesPerDir.rightward;
        rightwardLapStopTimes=representativeRightwardUnitInfo.lapStopTimesPerDir.rightward;
    end
    currSessionPosInfo=load(currPosInfoFilePath);
    
    %missing seq files?
    if(exist(currSeqDataFilePath,'file'))
        currThetaSeqData=load(currSeqDataFilePath);
    else
        continue
    end
    
    %numUnits=size(currThetaSeqData.avgSpikeTimesByCycle,2);
     if(isfield(currThetaSeqData,'totalNumFieldsInSession'))
         numFields=currThetaSeqData.totalNumFieldsInSession;
     else
         continue
     end
     numCycles=currThetaSeqData.numCycles;
     %sortedMeanCenterIdxPerDi=currThetaSeqData.sortedMeanCenterIdxPerDi;
     timeOfCycleTrough=currThetaSeqData.timeOfCycleTrough;
     peakSpikePhasePerCyclePerField=currThetaSeqData.peakSpikePhasePerCyclePerField;
     avgSpikeTimesPerCyclePerField=currThetaSeqData.avgSpikeTimesPerCyclePerField;
     
     spikeDataFieldIdxesPerField=currThetaSeqData.allSpikeDataFieldIdxesThisSession;
     
     allSpikeTimesPerCyclePerField=currThetaSeqData.allSpikeTimesByCyclePerField;
     
     allFieldIDsThisSession=currThetaSeqData.allFieldIDsThisSession;
     
     fieldCenterPerField=currThetaSeqData.fieldCenterPerField;
     [rightwardSortedFieldCenters, rightwardSortedFieldCenterIdxes]=sort(fieldCenterPerField);
     [leftwardSortedFieldCenters, leftwardSortedFieldCenterIdxes]=sort(fieldCenterPerField,'descend');

          
          fieldEndPerField=currThetaSeqData.fieldEndPerField;
          fieldStartPerField=currThetaSeqData.fieldStartPerField;
           
          posTimeAxis=currSessionPosInfo.positionTimeAxisSec;
          posPerTime=currSessionPosInfo.posAlongTrackPerTimeM;
          speedPerTime=abs(currSessionPosInfo.speedAlongTrackPerTimeMSec);
          maxTrackLength=currSessionPosInfo.trackLengthM;
          
          approxLapTime=10; %sec
          
          waveformsPerFieldPerSpike=cell(max(spikeDataFieldIdxesPerField),1);
          
           allTimeFracsPerSpikePerField=cell(max(spikeDataFieldIdxesPerField),1);
           firingRatePerPosPerField=cell(max(spikeDataFieldIdxesPerField),1);
           
          for sdi=1:length(spikeDataFieldIdxesPerField)
              currSpikeDataFieldIdx=spikeDataFieldIdxesPerField(sdi);
              %waveformsPerFieldPerSpike{currSpikeDataFieldIdx}=spikeDataPerField.thetaWaveformLFPzPerSpikePerField{currSpikeDataFieldIdx};
              if(useOmarAsymAdjustment)
                  waveformsPerFieldPerSpike{currSpikeDataFieldIdx}=spikeDataPerField.omarAsymAdjustedThetaWaveformLFPzPerSpikePerField{currSpikeDataFieldIdx};
              else
                  waveformsPerFieldPerSpike{currSpikeDataFieldIdx}=spikeDataPerField.adjustedThetaWaveformLFPzPerSpikePerField{currSpikeDataFieldIdx};
              end
              
              allTimeFracsPerSpikePerField{currSpikeDataFieldIdx}=spikeDataPerField.spikeTimeFracInFieldPerField{currSpikeDataFieldIdx};

             firingRatePerPosPerField{currSpikeDataFieldIdx}=spikeDataPerField.firingRatePerPosPerField{currSpikeDataFieldIdx};

          end
          spikeTimesInExpPerField=spikeDataPerField.spikeTimesInExpPerField;
          
          unitInfoPathPerField=spikeDataPerField.unitInfoPathPerField;
          
          if(~isempty(representativeLeftwardUnitInfo))
              representativeSessUnitInfo=representativeLeftwardUnitInfo;
          end
          
          if(~isempty(representativeRightwardUnitInfo))
              representativeSessUnitInfo=representativeRightwardUnitInfo;
          end
          
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %use adjusted cycle info
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          currSessName
          adjustedCycleInfoFilePath=getRegexFilePath(adjustedCycleDataDir,sprintf('*_%s_adjCycleLFPstruct.mat', currSessName));

          cycleInfoStruct=load(adjustedCycleInfoFilePath);
          refChStr=lower(getSubstrBtwnSubstrs(representativeSessUnitInfo.refChCycleInfoPath,sprintf('%s_',currSessName),'_6-12_minmax.mat'));

          cycleInfo=cycleInfoStruct.adjustedCycleInfoPerPyrCh.(refChStr);
          
          
          %cycleInfo=load(representativeSessUnitInfo.refChCycleInfoPath);
          
          if(useOmarAsymAdjustment)
              cycleInfoMinTimes=cycleInfo.fastMidTimes;
              cycleInfoMaxTimes=cycleInfo.fastStartTimes;
          else
              cycleInfoMinTimes=cycleInfo.cycleMinTimes;
              cycleInfoMaxTimes=cycleInfo.closestFastMaxTimes;
          end
          
          %approxLapTime=20; %sec
          
          approxTravTime=1;
    for ci=1:numCycles
        
        currCycleTime=timeOfCycleTrough(ci);

        [~,closestMaxIdx]=min(abs(cycleInfoMaxTimes-currCycleTime));
        if(cycleInfoMaxTimes(closestMaxIdx)>currCycleTime)
            closestMaxIdx=max(1,closestMaxIdx-1);
        end
        prevMaxIdx=closestMaxIdx;
        nextMaxIdx=min(length(cycleInfoMaxTimes),closestMaxIdx+1);
       
        
        currCycleStartTime=cycleInfoMaxTimes(prevMaxIdx);
        currCycleEndTime=cycleInfoMaxTimes(nextMaxIdx);
        
        currCycleDuration=currCycleEndTime-currCycleStartTime;
         
        currCyclePeakPhasePerField=peakSpikePhasePerCyclePerField(ci,:);
        currCycleSpikeTimesPerField=allSpikeTimesPerCyclePerField(ci,:);
        
        %currCycleMeanPhasePerField
        
        currCycleAllSpikeTimesPerField=allSpikeTimesPerCyclePerField(ci,:);
        
        %currCycleTimeFracPerField=
        
        %at least n fields in sequence
        if(sum(~isnan(currCyclePeakPhasePerField))<minNumFieldsInSeq)
            continue
        end
        
        
        currCycleMeanPhasePerField=NaN(size(currCyclePeakPhasePerField));
        for fci=1:length(currCyclePeakPhasePerField)
            currCycleFieldAllSpikePhases=[360*currCycleSpikeTimesPerField{fci}]/currCycleDuration;
            currCycleMeanPhasePerField(fci)=circMeanDeg(currCycleFieldAllSpikePhases);
         end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %get waveform of this cycle
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fieldIdxActiveThisCycle=find(~isnan(currCyclePeakPhasePerField));
        
       currCycleTimeFracPerSpikeDataFieldIdx=NaN(max(spikeDataFieldIdxesPerField),1);
        
        for w=1:length(fieldIdxActiveThisCycle)
            spikeDataFieldIdx=spikeDataFieldIdxesPerField(fieldIdxActiveThisCycle(w));
            waveformsPerSpike=waveformsPerFieldPerSpike{spikeDataFieldIdx};
            currFieldTimeFracsPerSpike=allTimeFracsPerSpikePerField{spikeDataFieldIdx};
            
            spikeTimesInExp=spikeTimesInExpPerField{spikeDataFieldIdx};
            
            currCycleSpikeTimesInExpIdxes=spikeTimesInExp>=currCycleStartTime & spikeTimesInExp<=currCycleEndTime;
            
            firstSpikeTimeInExpIdx=find(currCycleSpikeTimesInExpIdxes);
            
            if(~isempty(firstSpikeTimeInExpIdx))
                break
            end
        end
        
        
        firstSpikeTimeInExpIdx=find(currCycleSpikeTimesInExpIdxes);
        
        if(isempty(firstSpikeTimeInExpIdx))
            continue
        end
        
        firstSpikeTimeInExpIdx=firstSpikeTimeInExpIdx(1);
        
        currCycleWaveform=waveformsPerSpike{firstSpikeTimeInExpIdx};
        
        
        %firingRatePerPosPerFieldInCurrSeq=cell(length(fieldIdxActiveThisCycle),1);
        for w=1:length(fieldIdxActiveThisCycle)
            spikeDataFieldIdx=spikeDataFieldIdxesPerField(fieldIdxActiveThisCycle(w));
           
            currFieldTimeFracsPerSpike=allTimeFracsPerSpikePerField{spikeDataFieldIdx};
            %firingRatePerPosPerFieldInCurrSeq{w}=firingRatePerPosPerField{spikeDataFieldIdx};
            
            spikeTimesInExp=spikeTimesInExpPerField{spikeDataFieldIdx};
            
            currCycleSpikeTimesInExpIdxes=spikeTimesInExp>=currCycleStartTime & spikeTimesInExp<=currCycleEndTime;
            
            currCycleTimeFracPerSpikeDataFieldIdx(spikeDataFieldIdx)=nanmean(currFieldTimeFracsPerSpike(currCycleSpikeTimesInExpIdxes)); %all spikes should be same time frac
        end
        
      
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %convert from spikeDataFieldIdx to this session sequence data idx
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        currCycleTimeFracPerField=currCycleTimeFracPerSpikeDataFieldIdx(spikeDataFieldIdxesPerField);
        
        currCyclePresentFieldIdxes=find(~isnan(currCyclePeakPhasePerField));
        
        currCyclePresentFieldIDs=allFieldIDsThisSession(currCyclePresentFieldIdxes);
        
        currCyclePresentFieldDirections=NaN(size(currCyclePresentFieldIDs));
        
        allFieldsRightward=1;
        allFieldsLeftward=1;
        for fi=1:length(currCyclePresentFieldIDs)
            currFieldID=currCyclePresentFieldIDs(fi);
            [ui,di,dfi]=recoverFieldInfoFromID(currFieldID);
            
            if(di==1)
                allFieldsLeftward=0; %one rightward field
            elseif(di==2)
                allFieldsRightward=0; %one leftward field
            end
        end

       assert(~(allFieldsRightward && allFieldsLeftward));
        
       currDirSortedFieldCenterIdxes=NaN;
        if(allFieldsRightward)
            currCyclePeakPhasePerFieldInSpaceOrder=currCyclePeakPhasePerField(rightwardSortedFieldCenterIdxes);
            currCycleMeanPhasePerFieldInSpaceOrder=currCycleMeanPhasePerField(rightwardSortedFieldCenterIdxes);
            
            currCycleTimeFracPerFieldInSpaceOrder=currCycleTimeFracPerField(rightwardSortedFieldCenterIdxes);

            currCycleFieldCenterPerFieldInSpaceOrder=fieldCenterPerField(rightwardSortedFieldCenterIdxes);
            ratCurrCyclePos=interp1(posTimeAxis,posPerTime,currCycleTime);
            
            lapStartTimes=rightwardLapStartTimes;
            lapStopTimes=rightwardLapStopTimes;
            
            currSeqDirStr='rightward';
            
            currDirSortedFieldCenterIdxes=rightwardSortedFieldCenterIdxes;

        elseif(allFieldsLeftward)
            currCyclePeakPhasePerFieldInSpaceOrder=currCyclePeakPhasePerField(leftwardSortedFieldCenterIdxes);
            currCycleMeanPhasePerFieldInSpaceOrder=currCycleMeanPhasePerField(leftwardSortedFieldCenterIdxes);
            
            currCycleTimeFracPerFieldInSpaceOrder=currCycleTimeFracPerField(leftwardSortedFieldCenterIdxes);
            
            currCycleFieldCenterPerFieldInSpaceOrder=fieldCenterPerField(leftwardSortedFieldCenterIdxes);
            
            ratCurrCyclePos=interp1(posTimeAxis,maxTrackLength-posPerTime,currCycleTime);
            
            lapStartTimes=leftwardLapStartTimes;
            lapStopTimes=leftwardLapStopTimes;
            
            currSeqDirStr='leftward';
            
            currDirSortedFieldCenterIdxes=leftwardSortedFieldCenterIdxes;
        end
        
        currCycleSpikeTimesPerFieldInSpaceOrder=currCycleSpikeTimesPerField(currDirSortedFieldCenterIdxes);
        
        fieldStartPerFieldInSpaceOrder=fieldStartPerField(currDirSortedFieldCenterIdxes);
        fieldEndPerFieldInSpaceOrder=fieldEndPerField(currDirSortedFieldCenterIdxes);
        
        currCycleActiveFieldOrderedIdxes=~isnan(currCyclePeakPhasePerFieldInSpaceOrder);
        
        currCyclePhaseSeq=currCyclePeakPhasePerFieldInSpaceOrder(currCycleActiveFieldOrderedIdxes);
        currCycleMeanPhaseSeq=currCycleMeanPhasePerFieldInSpaceOrder(currCycleActiveFieldOrderedIdxes);
        
        
        currCycleTimeFracSeq=currCycleTimeFracPerFieldInSpaceOrder(currCycleActiveFieldOrderedIdxes);
        

          spikeDataFieldIdxesPerFieldInSpaceOrder=spikeDataFieldIdxesPerField(currDirSortedFieldCenterIdxes);
        
        
        if(allFieldsLeftward)
            currCycleFieldCenterSeq=maxTrackLength-currCycleFieldCenterPerFieldInSpaceOrder(currCycleActiveFieldOrderedIdxes); %1-distances for leftward runs
            currCycleFieldDistStartSeq=maxTrackLength-fieldStartPerFieldInSpaceOrder(currCycleActiveFieldOrderedIdxes);
            currCycleFieldDistEndSeq=maxTrackLength-fieldEndPerFieldInSpaceOrder(currCycleActiveFieldOrderedIdxes);
            
            lapsStartTimes=leftwardLapStartTimes;
            lapStopTimes=leftwardLapStopTimes;
           
        else
             currCycleFieldCenterSeq=currCycleFieldCenterPerFieldInSpaceOrder(currCycleActiveFieldOrderedIdxes);
            currCycleFieldDistStartSeq=fieldStartPerFieldInSpaceOrder(currCycleActiveFieldOrderedIdxes);
            currCycleFieldDistEndSeq=fieldEndPerFieldInSpaceOrder(currCycleActiveFieldOrderedIdxes);
            
            lapsStartTimes=rightwardLapStartTimes;
            lapStopTimes=rightwardLapStopTimes;
        end
        
        
        [~,closestsLapStartIdx]=min(abs(currCycleTime-lapsStartTimes));
        if(currCycleTime<lapsStartTimes(closestsLapStartIdx))
            closestsLapStartIdx=closestsLapStartIdx+1;
        end
        
        currLapIdx=min(closestsLapStartIdx,length(lapsStartTimes));
        
        currCycleFieldWidthSeq=currCycleFieldDistEndSeq-currCycleFieldDistStartSeq;
        avgFieldWidth=nanmean(currCycleFieldWidthSeq);
        %currCycleFieldCenterSeqAroundCenter=currCycleFieldCenterSeq-nanmean(currCycleFieldCenterSeq);
        
        thisCycleSpikeDataFieldIdxes=spikeDataFieldIdxesPerFieldInSpaceOrder(currCycleActiveFieldOrderedIdxes);
        
        currThetaSeqFieldIdxStr='fieldSeq';
        
        for ai=1:length(thisCycleSpikeDataFieldIdxes)
            currFieldIdx=thisCycleSpikeDataFieldIdxes(ai);
            currThetaSeqFieldIdxStr=[currThetaSeqFieldIdxStr '_' num2str(currFieldIdx)];
        end
        
        middleFieldSpikeDataIdx=thisCycleSpikeDataFieldIdxes(floor(length(thisCycleSpikeDataFieldIdxes)/2));
        
        currCycleSeqAllSpikeTimesPerField=currCycleSpikeTimesPerFieldInSpaceOrder(currCycleActiveFieldOrderedIdxes);
        
        currCycleSeqAllSpikePhasesPerField=cell(size(currCycleSeqAllSpikeTimesPerField));
        for ai=1:length(currCycleSeqAllSpikeTimesPerField)
            currCycleSeqAllSpikePhasesPerField{ai}=360*(currCycleSeqAllSpikeTimesPerField{ai}/currCycleDuration);
        end
        
        %{
        fieldCountInSeq=0;
        currCycleAllPhaseSeq=[];
        for ai=1:length(currCycleActiveFieldOrderedIdxes)
            if(currCycleActiveFieldOrderedIdxes(ai))
                fieldCountInSeq=fieldCountInSeq+1;
                currCycleAllPhaseSeq{fieldCountInSeq}=currCycleSeqAllSpikePhasesPerField{ai};
            end
        
        end
        %}
        
        middleFieldThetaData=spikeDataPerField.thetaDataPerField{middleFieldSpikeDataIdx};
        
        middleFieldMeanSpeedsAcrossLaps=middleFieldThetaData.allLapInFieldMeanSpeedPerTraversal;
        
        allSpikePhasesPerField=cell(length(currCycleSeqAllSpikeTimesPerField),1);
        allSpikeTimeFracsInFieldPerField=cell(length(currCycleSeqAllSpikeTimesPerField),1);
        allSpikePosPerField=cell(length(currCycleSeqAllSpikeTimesPerField),1);
        for ai=1:length(currCycleSeqAllSpikeTimesPerField)
            currFieldSpikeDataIdx=thisCycleSpikeDataFieldIdxes(ai);
            
            allSpikePhasesPerField{ai}=spikeDataPerField.spikePhasesPerField{currFieldSpikeDataIdx};
            allSpikeTimesThisField=spikeDataPerField.spikeTimesInExpPerField{currFieldSpikeDataIdx};
            
            allSpikeTimeFracsInFieldPerField{ai}=spikeDataPerField.spikeTimeFracInFieldPerField{currFieldSpikeDataIdx};
            
            
            if(allFieldsRightward)
                allSpikePosPerField{ai}=interp1(posTimeAxis,posPerTime,allSpikeTimesThisField);
            elseif(allFieldsLeftward)
                allSpikePosPerField{ai}=interp1(posTimeAxis,maxTrackLength-posPerTime,allSpikeTimesThisField);

            end
        end
        
        %currMeanTravSpeed=middleFieldMeanSpeedsAcrossLaps(currLapIdx);
        
         %approxCurrLapStartTime=currCycleTime-approxLapTime/2;
         %approxCurrLapEndTime=currCycleTime+approxLapTime/2;
         
         currSeqStartPos=min(currCycleFieldDistStartSeq);
         currSeqEndPos=max(currCycleFieldDistEndSeq);
         
         currLapStartTime=lapStartTimes(currLapIdx);
         currLapStopTime=lapStopTimes(currLapIdx);
         
        
         currLapPosTimeIdxes=posTimeAxis>=currLapStartTime & posTimeAxis<=currLapStopTime;
         
         %currLapPosPerTime=posPerTime(currLapPosTimeIdxes);
         
         if(allFieldsLeftward)
             standardizedPosPerTime=maxTrackLength-posPerTime;
         elseif(allFieldsRightward)
             standardizedPosPerTime=posPerTime;
         end
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %get fraction of field traversed per field at this cycle
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %ratCurrCyclePos
         currCycleDistFracSeq=NaN(size(currCycleFieldDistStartSeq));
         for ai=1:length(currCycleFieldDistStartSeq)
             
            currFieldStartDist=currCycleFieldDistStartSeq(ai);
            currFieldEndDist=currCycleFieldDistEndSeq(ai);
            currFieldWidth=currCycleFieldWidthSeq(ai);
            
            currFieldDistFrac=(ratCurrCyclePos-currFieldStartDist)/(currFieldWidth);
            
            currCycleDistFracSeq(ai)=currFieldDistFrac; 
         end
         
         disp('')
         
         
         %{
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %get start and stop times per field for this traversal
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
         for ai=1:length(currCycleFieldStartSeq)
             %[startBoundTime,endBoundTime,withinBoundsTimeIdxes] = ...
                 %getTimeBoundsOutwardFromSeedAndPosBounds(lowestPossibleTime,highestPossibleTime,timeAxis,posPerTime,centralSeedTime,lowerPosBound,upperPosBound)

                currFieldPosStart= currCycleFieldStartSeq(ai);
                currFieldPosEnd= currCycleFieldEndSeq(ai);
                
             [currTravStartTimeCurrField,endBoundTime,withinBoundsTimeIdxes] = ...
                 getTimeBoundsOutwardFromSeedAndPosBounds(lowestPossibleTime,highestPossibleTime,posTimeAxis,standardizedPosPerTime,currCycleTime,currFieldPosStart,currFieldPosEnd)

             currTravStartTimePerField{ai}=currTravStartTimeCurrField;
             
         end
         %}
             
             
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %step backward from current cycle seed to find start of trav
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         [~,currLapPosTimeIdxStart]=min(abs(currLapStartTime - posTimeAxis));
         [~,currTravPosTimeIdxStart]=min(abs(currCycleTime - posTimeAxis));
         while(standardizedPosPerTime(currTravPosTimeIdxStart)>currSeqStartPos)
             currTravPosTimeIdxStart=currTravPosTimeIdxStart-1;
             
             if(currTravPosTimeIdxStart<currLapPosTimeIdxStart)
                 break
             end
         end
         
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %step forward from current cycle seed to find end of trav
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         [~,currLapPosTimeIdxEnd]=min(abs(currLapStopTime - posTimeAxis));
         [~,currTravPosTimeIdxEnd]=min(abs(currCycleTime - posTimeAxis));
         while(standardizedPosPerTime(currTravPosTimeIdxEnd)<currSeqEndPos)
             currTravPosTimeIdxEnd=currTravPosTimeIdxEnd+1;
             
             if(currTravPosTimeIdxEnd>currLapPosTimeIdxEnd)
                 break
             end
         end
         
         
         %currTravPosTimeIdxes=currLapPosPerTime>=currSeqStartPos & currLapPosPerTime<=currSeqEndPos;

         currTravPosTimeIdxes=currTravPosTimeIdxStart:currTravPosTimeIdxEnd;
         
         currTravStartTimeInExp=posTimeAxis(currTravPosTimeIdxStart);
         currTravEndTimeInExp=posTimeAxis(currTravPosTimeIdxEnd);
         currTimeInExpWithinTrav=posTimeAxis(currTravPosTimeIdxes);
         
         currTravStartDist=standardizedPosPerTime(currTravPosTimeIdxStart);
         currTravEndDist=standardizedPosPerTime(currTravPosTimeIdxEnd);
         currTravDistPerTime=standardizedPosPerTime(currTravPosTimeIdxes);
         
         %currTravPosPerTime=currLapPosPerTime(currTravPosTimeIdxes);
                    
        currMeanTravSpeed=nanmean(speedPerTime(currTravPosTimeIdxes));
        currTravSpeedTrace=speedPerTime(currTravPosTimeIdxes);
        
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %get time within field traversal for this cycle per field
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         startTimeInFieldPerField=NaN(size(currCycleFieldDistStartSeq));
         endTimeInFieldPerField=NaN(size(currCycleFieldDistStartSeq));
         centerTimeOfFieldPerField=NaN(size(currCycleFieldDistStartSeq));
         
         for fti=1:length(currCycleFieldDistStartSeq)
             currFieldStartPos=currCycleFieldDistStartSeq(fti);
             currFieldEndPos=currCycleFieldDistEndSeq(fti);
             
             currFieldCenterPos=currCycleFieldCenterSeq(fti);
             
             %interpIdxes=~isnan(currTravDistPerTime) & ~isnan(currTimeInExpWithinTrav);
             %currFieldStartTime=interp1(currTravDistPerTime(interpIdxes),currTimeInExpWithinTrav(interpIdxes),currFieldStartPos);
            filledTimeInTrav=fillmissing(currTravDistPerTime,'linear');
            filledDistInTrav=fillmissing(currTimeInExpWithinTrav,'linear');
            
            %take only first of repeat value indexes
            [~,uniqIdxes]=unique(filledTimeInTrav);
             
             currFieldStartTime=interp1(filledTimeInTrav(uniqIdxes),filledDistInTrav(uniqIdxes),currFieldStartPos,'linear','extrap');
             currFieldEndTime=interp1(filledTimeInTrav(uniqIdxes),filledDistInTrav(uniqIdxes),currFieldEndPos,'linear','extrap');

             
             currFieldCenterTime=interp1(filledTimeInTrav(uniqIdxes),filledDistInTrav(uniqIdxes),currFieldCenterPos,'linear','extrap');
      
             
             startTimeInFieldPerField(fti)=currFieldStartTime;
             endTimeInFieldPerField(fti)=currFieldEndTime;
             centerTimeOfFieldPerField(fti)=currFieldCenterTime;
         end
         
        currCycleTimeSinceFieldStartPerFieldInSec=currCycleTime-startTimeInFieldPerField; %adjust each by field center time....
        
        currCycleTimeToFieldEndPerFieldInSec=endTimeInFieldPerField-currCycleTime;
        
        %keep track of time since field start (a la Hopfield)
        %currCycleTimeOfFieldInSeqPerFieldInSec=(startTimeInFieldPerField-currTravStartTimeInExp);
        currCycleTimeOfFieldInSeqPerFieldInSec=(centerTimeOfFieldPerField-currTravStartTimeInExp);
        
       %time since start of field + time offset of this field within sequence
        %currCycleTimeSinceStartOfFieldPlusFieldTimeOffestPerField=currCycleTimeSinceFieldStartPerFieldInSec+currCycleTimeOfFieldInSeqPerFieldInSec;
       
        %time within field relative to sequence interval
        
        %currCycleTimeInFieldRelativeToSeqPerFieldInSec=(currCycleTime-startTimeInFieldPerField)-currTravStartTimeInExp;
        
        currCycleTimeSinceTravStartPerFieldInSec=(currCycleTime-currTravStartTimeInExp);
        
        currCycleTotalTraversalDurationSec=currTravEndTimeInExp-currTravStartTimeInExp;
        
        %currCycleTimeToTypicalFieldEndPerFieldInSec=typicalFieldDurationPerField-currCycleTime;
        
        currTravFieldDurationPerFieldInSec=endTimeInFieldPerField-startTimeInFieldPerField;
        

        %angDiffPerField=angdiffDeg(currCyclePhaseSeq);
        
        
        if(onlyCorrectOrder)
            
            %if not increasing, on average
            %if(mean(angDiffPerField)<0)
            %if not always increasing
            %{
            if(min(angDiffPerField)<0)
                continue
            end
            %}
            
             if(min(diff(currCyclePhaseSeq))<0)
                continue
            end
            
            %[rho,p,~, ~]= getCircCorrCoeff(currCycleFieldCenterSeq,currCyclePhaseSeq);
        end
        
        if(onlySimilarSizeFields)
            maxFieldWidth=max(currCycleFieldWidthSeq);
            minFieldWidth=min(currCycleFieldWidthSeq);
            
            if(minFieldWidth/maxFieldWidth<minAllowableFieldRatio)
                continue
            end
        end
        
        approxTrajectoryLength=range(currCycleFieldCenterSeq);
        meanFieldWidth=nanmean(currCycleFieldWidthSeq);
        
        trajLengthFieldWidthRatio=approxTrajectoryLength/meanFieldWidth;
        
        minStartDistFieldWidthRatio=diff(currCycleFieldDistStartSeq)/meanFieldWidth;
        
        minEndDistFieldWidthRatio= diff(currCycleFieldDistEndSeq)/meanFieldWidth;
        
        if(enforceMinTrajectoryLength)
            if(trajLengthFieldWidthRatio<minTrajLengthWidthRatio)
                continue
            end 
            
            %control start and end separation 
            %{
            if(minStartDistFieldWidthRatio<minTrajLengthWidthRatio)
                continue
            end 
            
            if(minEndDistFieldWidthRatio<minTrajLengthWidthRatio)
                continue
            end 
            %}
        end
        
        
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %initialize sequence info struct for this field sequence or
        %increment occurence counter if same as previous
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if(~isfield(thetaSequenceInfoPerSequence,currThetaSeqFieldIdxStr))
           
            thetaSequenceInfoPerSequence.(currThetaSeqFieldIdxStr)=[];
             thetaSequenceInfoPerSequence.(currThetaSeqFieldIdxStr).numOccurences=1;
            thetaSequenceInfoPerSequence.(currThetaSeqFieldIdxStr).seqDirStr=currSeqDirStr;
            thetaSequenceInfoPerSequence.(currThetaSeqFieldIdxStr).avgFieldWidth=avgFieldWidth;
             
            thetaSequenceInfoPerSequence.(currThetaSeqFieldIdxStr).spikeDataFieldIdxList=thisCycleSpikeDataFieldIdxes;
            thetaSequenceInfoPerSequence.(currThetaSeqFieldIdxStr).fieldCenterSeq=currCycleFieldCenterSeq;
            thetaSequenceInfoPerSequence.(currThetaSeqFieldIdxStr).lapStartTimesUsedForLapIdx=lapStartTimes;
            thetaSequenceInfoPerSequence.(currThetaSeqFieldIdxStr).lapStopTimesUsedForLapIdx=lapStopTimes;
            
            thetaSequenceInfoPerSequence.(currThetaSeqFieldIdxStr).lapStopTimesUsedForLapIdx=lapStopTimes;
            thetaSequenceInfoPerSequence.(currThetaSeqFieldIdxStr).lapStopTimesUsedForLapIdx=lapStopTimes;

            thetaSequenceInfoPerSequence.(currThetaSeqFieldIdxStr).meanPhaseSeqPerOcc=[];
            thetaSequenceInfoPerSequence.(currThetaSeqFieldIdxStr).allPhaseSeqPerOcc=[];
            thetaSequenceInfoPerSequence.(currThetaSeqFieldIdxStr).timeSinceFieldStartSeqPerOcc=[];
            thetaSequenceInfoPerSequence.(currThetaSeqFieldIdxStr).timeToFieldEndSeqPerOcc=[];
            thetaSequenceInfoPerSequence.(currThetaSeqFieldIdxStr).timeOfFieldInSeqPerOccPerFieldInSec=[];
            thetaSequenceInfoPerSequence.(currThetaSeqFieldIdxStr).totalTraversalDurationPerOcc=[];
            thetaSequenceInfoPerSequence.(currThetaSeqFieldIdxStr).timeSinceTravStartEveryCyclePerOcc=[];

            
            thetaSequenceInfoPerSequence.(currThetaSeqFieldIdxStr).fieldDurationSeqPerOcc=[];
            
            
            thetaSequenceInfoPerSequence.(currThetaSeqFieldIdxStr).peakPhaseSeqPerOcc=[];
            thetaSequenceInfoPerSequence.(currThetaSeqFieldIdxStr).fieldDistFracSeqPerOcc=[];
            thetaSequenceInfoPerSequence.(currThetaSeqFieldIdxStr).fieldTimeFracSeqPerOcc=[];
            thetaSequenceInfoPerSequence.(currThetaSeqFieldIdxStr).cycleNumPerOcc=[];
            thetaSequenceInfoPerSequence.(currThetaSeqFieldIdxStr).positionPerOcc=[];
            thetaSequenceInfoPerSequence.(currThetaSeqFieldIdxStr).lapNumPerOcc=[];
            thetaSequenceInfoPerSequence.(currThetaSeqFieldIdxStr).cycleTimeInExpPerOcc=[];
            thetaSequenceInfoPerSequence.(currThetaSeqFieldIdxStr).travSpeedTracePerOcc=[];
            thetaSequenceInfoPerSequence.(currThetaSeqFieldIdxStr).travAvgSpeedPerOcc=[];
            thetaSequenceInfoPerSequence.(currThetaSeqFieldIdxStr).thetaCycleWaveform=[];
        else
            thetaSequenceInfoPerSequence.(currThetaSeqFieldIdxStr).numOccurences=thetaSequenceInfoPerSequence.(currThetaSeqFieldIdxStr).numOccurences+1;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %add sequence info struct for this sequence occurence.
        %info desired per sequence:
        %-spike data field idxes
        %-field center sequence
        %-phase sequence per occurence
        %-cycle number per occurence
        %-position per occurence
        %-time in experiment per occurence
        %-avg speed per occurence
        %-speed trace per occurence
        %-lap number per occurence
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         oi=thetaSequenceInfoPerSequence.(currThetaSeqFieldIdxStr).numOccurences;
         
         thetaSequenceInfoPerSequence.(currThetaSeqFieldIdxStr).fieldTimeFracSeqPerOcc{oi}=currCycleTimeFracSeq;
         thetaSequenceInfoPerSequence.(currThetaSeqFieldIdxStr).fieldDistFracSeqPerOcc{oi}=currCycleDistFracSeq;
         
           thetaSequenceInfoPerSequence.(currThetaSeqFieldIdxStr).meanPhaseSeqPerOcc{oi}=currCycleMeanPhaseSeq;
           thetaSequenceInfoPerSequence.(currThetaSeqFieldIdxStr).allPhaseSeqPerOcc{oi}=currCycleSeqAllSpikePhasesPerField;
            thetaSequenceInfoPerSequence.(currThetaSeqFieldIdxStr).timeSinceFieldStartSeqPerOcc{oi}=currCycleTimeSinceFieldStartPerFieldInSec;
             thetaSequenceInfoPerSequence.(currThetaSeqFieldIdxStr).timeToFieldEndSeqPerOcc{oi}=currCycleTimeToFieldEndPerFieldInSec;
             thetaSequenceInfoPerSequence.(currThetaSeqFieldIdxStr).timeSinceTravStartEveryCyclePerOcc{oi}=currCycleTimeSinceTravStartPerFieldInSec;
            thetaSequenceInfoPerSequence.(currThetaSeqFieldIdxStr).fieldDurationSeqPerOcc{oi}=currTravFieldDurationPerFieldInSec;
            thetaSequenceInfoPerSequence.(currThetaSeqFieldIdxStr).totalTraversalDurationPerOcc{oi}=currCycleTotalTraversalDurationSec;
           
            thetaSequenceInfoPerSequence.(currThetaSeqFieldIdxStr).timeOfFieldInSeqPerOccPerFieldInSec{oi}=currCycleTimeOfFieldInSeqPerFieldInSec;
         
        thetaSequenceInfoPerSequence.(currThetaSeqFieldIdxStr).peakPhaseSeqPerOcc{oi}=currCyclePhaseSeq;
        thetaSequenceInfoPerSequence.(currThetaSeqFieldIdxStr).cycleNumPerOcc{oi}=ci;
            thetaSequenceInfoPerSequence.(currThetaSeqFieldIdxStr).positionPerOcc{oi}=ratCurrCyclePos;
            thetaSequenceInfoPerSequence.(currThetaSeqFieldIdxStr).lapNumPerOcc{oi}=currLapIdx;
            thetaSequenceInfoPerSequence.(currThetaSeqFieldIdxStr).cycleTimeInExpPerOcc{oi}=currCycleTime;
            thetaSequenceInfoPerSequence.(currThetaSeqFieldIdxStr).travSpeedTracePerOcc{oi}=currTravSpeedTrace;
            thetaSequenceInfoPerSequence.(currThetaSeqFieldIdxStr).travAvgSpeedPerOcc{oi}=currMeanTravSpeed;
            thetaSequenceInfoPerSequence.(currThetaSeqFieldIdxStr).thetaCycleWaveform{oi}=currCycleWaveform;
            
            firingRatePerPosPerFieldInSeq=[];
            for ai=1:length(currCycleFieldCenterSeq)
              
                spikeDataFieldIdx= thisCycleSpikeDataFieldIdxes(ai);
                firingRatePerPos=firingRatePerPosPerField{spikeDataFieldIdx};
                
        
                
                if(allFieldsLeftward)
                     distBins=maxTrackLength - firingRatePerPos(:,1);
                else
                     distBins=firingRatePerPos(:,1);
                end
                
                distBins(distBins>currTravEndDist | distBins<currTravStartDist)=NaN; 
                
                outIdxes=distBins>currTravEndDist | distBins<currTravStartDist;
                
                firingRatePerPos(outIdxes,1)=NaN;
                
                firingRatePerPosPerFieldInSeq{ai}=firingRatePerPos;

            end

           thetaSequenceInfoPerSequence.(currThetaSeqFieldIdxStr).firingRatePerPosPerFieldInSeq=firingRatePerPosPerFieldInSeq;
        
       
                if(showPlots)
                    
             fH=figure;
             
            %plot(currCyclePhaseSeq,ones(size(currCyclePhaseSeq)),'k.')
            currCycleFieldCount=length(currCyclePhaseSeq);
            
            subplot(2,2,1)
           distColors=jet(length(currCycleFieldCenterSeq))
           for ai=1:length(currCycleFieldCenterSeq)
               %plot(currCycleFieldCenterSeq(ai),currCyclePhaseSeq(ai),'.','Color',distColors(ai,:),'MarkerSize',50)
               %plot(currCyclePhaseSeq(ai),currCycleFieldCenterSeq(ai),'.','Color',distColors(ai,:),'MarkerSize',50)
               %plot(currCyclePhaseSeq(ai),currCycleTimeFracSeq(ai),'.','Color',distColors(ai,:),'MarkerSize',50)
               plot(currCyclePhaseSeq(ai),1-currCycleDistFracSeq(ai),'.','Color',distColors(ai,:),'MarkerSize',50)

               
               hold on
               %plot([currCycleFieldStartSeq(ai) currCycleFieldEndSeq(ai)],[currCyclePhaseSeq(ai) currCyclePhaseSeq(ai)],'-','Color',distColors(ai,:),'LineWidth',5)
               %plot([currCyclePhaseSeq(ai) currCyclePhaseSeq(ai)],[currCycleFieldStartSeq(ai) currCycleFieldEndSeq(ai)],'-','Color',distColors(ai,:),'LineWidth',5)
 
               %vline(
           end
           colormap(jet)
           cb=colorbar
            ylabel(cb,'Behavioral field order')
            caxis([1 length(currCycleFieldCenterSeq)])

            %hold on
            %xlim([0 maxTrackLength])
            distRange=range([min(currCycleFieldDistStartSeq) max(currCycleFieldDistEndSeq)]);
            
            %{
            xlim([min(currCycleFieldStartSeq)-distRange*0.2 max(currCycleFieldEndSeq)+distRange*0.2])
            xlabel('Distance along track (m)')
            %ylim([-0.5 0.5])
            ylim([0 360])
            ylabel('Theta phase (deg)')
            %}
            
            %ylim([min(currCycleFieldStartSeq)-distRange*0.1 max(currCycleFieldEndSeq)+distRange*0.1])
            %ylabel('Distance along track (m)')
            
             %ylim([min(currCycleTimeFracSeq)-timeRange*0.1 max(currCycleTimeFracSeq)+distRange*0.1])
             ylim([-0.2 1.2])
            %ylabel('Time in respective field (frac)')
            ylabel('Distance to end of respective field (frac)')
            
            %ylim([-0.5 0.5])
            xlim([0 360])
            xlabel('Theta phase (deg)')
            
             
            %pL=vline(ratCurrCyclePos,'k--',3)
          
            %pL=hline(ratCurrCyclePos,'k--',3)
            hold off
            
            title(sprintf('%s, single theta sequence, cycle %d',currSessName, ci))
            %legend([pL pF],{'current position','active field'},'Location','northwest')%'eastoutside')
            %legend([pL],{'curr. position'},'Location','northwest')%'eastoutside')
            %legend boxoff
            axis square
            box off
            
            subplot(2,2,4)
            %plot(posPerTime(currLapPosTimeIdxes),speedPerTime(currLapPosTimeIdxes),'b-')
            
            timeVals=posTimeAxis(currTravPosTimeIdxes)-currCycleTime;
            if(allFieldsLeftward)
                %plotWithColorVec(maxTrackLength - posPerTime(currLapPosTimeIdxes),speedPerTime(currLapPosTimeIdxes),timeVals,'jet',5)
                %plotWithColorVec(speedPerTime(currLapPosTimeIdxes),maxTrackLength - posPerTime(currLapPosTimeIdxes),timeVals,'jet',5)
                %plotWithColorVec(speedPerTime(currTravPosTimeIdxes),maxTrackLength - posPerTime(currTravPosTimeIdxes),timeVals,'jet',5)
                %plotWithColorVec(maxTrackLength - posPerTime(currTravPosTimeIdxes),speedPerTime(currTravPosTimeIdxes),timeVals,'jet',5)
                plot(maxTrackLength - posPerTime(currTravPosTimeIdxes),speedPerTime(currTravPosTimeIdxes),'k-','LineWidth',3)

            else
               %plotWithColorVec(posPerTime(currLapPosTimeIdxes),speedPerTime(currLapPosTimeIdxes),timeVals,'jet',5)
                %plotWithColorVec(speedPerTime(currTravPosTimeIdxes),posPerTime(currTravPosTimeIdxes),timeVals,'jet',5)
                %plotWithColorVec(posPerTime(currTravPosTimeIdxes),speedPerTime(currTravPosTimeIdxes),timeVals,'jet',5)
                 plot(posPerTime(currTravPosTimeIdxes),speedPerTime(currTravPosTimeIdxes),'k-','LineWidth',3)

            end
            hold on
            
            for ai=1:length(currCycleFieldCenterSeq)
              
                spikeDataFieldIdx= thisCycleSpikeDataFieldIdxes(ai);
                firingRatePerPos=firingRatePerPosPerField{spikeDataFieldIdx};
                
                yyaxis right
                
                if(allFieldsLeftward)
                     distBins=maxTrackLength - firingRatePerPos(:,1);
                else
                     distBins=firingRatePerPos(:,1);
                end
                
                distBins(distBins>currTravEndDist | distBins<currTravStartDist)=NaN; 
                
                outIdxes=distBins>currTravEndDist | distBins<currTravStartDist;
                
                firingRatePerPos(outIdxes,1)=NaN;
                
               plot(distBins,firingRatePerPos(:,2),'-','Color',distColors(ai,:),'LineWidth',3)

               
               hold on

            end
            ylabel('Spatial firing rate (Hz)')
             yyaxis left
             hold on
             
           thetaSequenceInfoPerSequence.(currThetaSeqFieldIdxStr).firingRatePerPosInTrav=firingRatePerPos;
            %{
            cb=colorbar
            ylabel(cb,'Time from curr. cycle in trav. (sec)')
            %}
            
                  colormap(jet)
           cb=colorbar
            ylabel(cb,'Behavioral field order')
            caxis([1 length(currCycleFieldCenterSeq)])
            
            
            %{
            xlim([min(currCycleFieldStartSeq)-distRange*0.2 max(currCycleFieldEndSeq)+distRange*0.2])
            xlabel('Distance along track (m)')
              ylim([0 1])
            ylabel('Running speed (m/s)')
            %}
            %{
             ylim([min(currCycleFieldStartSeq)-distRange*0.1 max(currCycleFieldEndSeq)+distRange*0.1])
             ylabel('Distance along track (m)')
              xlim([0 maxRunSpeed])
            xlabel('Running speed in curr. traversal (m/s)')
            %}
              xlim([min(currCycleFieldDistStartSeq)-distRange*0.1 max(currCycleFieldDistEndSeq)+distRange*0.1])
             xlabel('Distance along track (m)')
              ylim([0 maxRunSpeed])
            ylabel('Running speed in curr. traversal (m/s)')
            vline(ratCurrCyclePos,'k--',3)
            hold on
            
           
            %pL=hline(ratCurrCyclePos,'k--',3);
            
            %legend(pL,'curr. position','Location','northeast')
            %legend boxoff
             %title(sprintf('running speed across traversal'))
             
            axis square
            hold off
            
            %axes('Position',[.78 .38 .075 .075])
            axes('Position',[.78 .58 .075 .075])
            speedEdges=linspace(0,maxRunSpeed,21);
            histogram(middleFieldMeanSpeedsAcrossLaps,speedEdges)
            hold on
           
            
             pS=vline(currMeanTravSpeed,'k--',2);
             
             %legend([pS],{'curr. speed'},'Location','northeast')%'eastoutside')
             %legend boxoff
            
             xlim([0 maxRunSpeed])
            xlabel('trav. speed (m/s)')
            
            ylabel('lap count')
             axis tight
             axis square
             box off
            
            
            subplot(2,2,3)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %waveform and local raster
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            phaseBins=linspace(0,360,length(currCycleWaveform(:,1)));
            plot(phaseBins,currCycleWaveform(:,2),'k-','LineWidth',1)
            hold on
            
            
            tickLineWidth=3;
            for ai=1:length(currCycleFieldCenterSeq)
                %[pH]=plotRasterStyle(eventTimes,rowNum,startTick,endTick,tickColor,tickLineWidth)
                %plotRasterStyle(currCycleSeqAllSpikePhasesPerField{ai},ai-3,NaN,NaN,distColors(ai,:),tickLineWidth)
                plotRasterStyle(currCycleSeqAllSpikePhasesPerField{ai},0,NaN,NaN,distColors(ai,:),tickLineWidth)

            end
            
            
            
            ylim([-2.5 2.5])
            xlim([0 360])
            xlabel('Theta phase (deg)')
            ylabel('raw LFP (Z)')
            
            colormap(jet)
           cb=colorbar
            ylabel(cb,'Behavioral field order')
            caxis([1 length(currCycleFieldCenterSeq)])
            
            axis square
            
            %subplot(2,2,2)
            %{
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %avg speed distribution across laps
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes('Position',[.7 .7 .2 .2])
            speedEdges=linspace(0,maxRunSpeed,21);
            histogram(middleFieldMeanSpeedsAcrossLaps,speedEdges)
            hold on
            
             pS=vline(currMeanTravSpeed,'k--',3);
             
             legend([pS],{'curr. speed'},'Location','northeast')%'eastoutside')
             legend boxoff
            
             xlim([0 maxRunSpeed])
            xlabel('Mean running speed in traversal (m/s)')
            
            ylabel('lap count')
            %}
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %bigger picture place fields
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %{
            for ai=1:length(currCycleFieldCenterSeq)
                
                %[pH]=plotRasterStyle(eventTimes,rowNum,startTick,endTick,tickColor,tickLineWidth)
                plot(allSpikePhasesPerField{ai},allSpikePosPerField{ai},'.','Color',distColors(ai,:),'MarkerSize',9)
                
                hold on
                %plot(allSpikePhasesPerField{ai},allSpikePosPerField{ai},'ko','MarkerSize',5)
            end
       
            
            %ylim([0 maxTrackLength])
           ylim([min(currCycleFieldStartSeq)-distRange*0.1 max(currCycleFieldEndSeq)+distRange*0.1])

            xlim([0 360])
            ylabel('Distance along track (m)')
            xlabel('Theta phase (deg)')
            
            pL=hline(ratCurrCyclePos,'k--',1);
            
            title('Phase precession across all laps')
            
             colormap(jet)
           cb=colorbar
            ylabel(cb,'Behavioral field order')
            caxis([1 length(currCycleFieldCenterSeq)])
            
            %legend([pL],{'curr. position'},'Location','northwest')%'eastoutside')

            %legend boxoff
            
            axis square
            %}
            
            %title(sprintf('running speed across laps'))
            
            setFigFontTo(16)
            %maxFigHalfHalfWidth
            %maxFigMukkaalWidth
            maxFig
            box off
            
            if(strcmp(currSessName,'ec013.626') && ci== 3837)
                
                disp('')
                
            end
            if(strcmp(currSessName,'ec013.626') && ci== 2277)
                
                disp('')
                
            end
            
            
            saveas(gcf,fullfile(saveThetaSeqImgDir,sprintf('%s_ThetaCycle%d_%dFields.png',currSessName, ci,currCycleFieldCount)))
            %saveEPS(fullfile(saveThetaSeqImgDir,sprintf('%s_ThetaCycle%d_%dFields',currSessName, ci,currCycleFieldCount)))
            %drawnow
            %pause(0.3)
            close all
          end
  

           
    end
    
  
    %thetaSequenceInfoPerSequenceWith4Fields=thetaSequenceInfoPerSequence;
    save(currSeqDataFilePath,'thetaSequenceInfoPerSequence','-append')
    %save(currSeqDataFilePath,'thetaSequenceInfoPerSequenceWith4Fields','-append')
end
        
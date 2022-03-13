close all;  clc
%clear all
clearvars -except spikeDataPerField
%dataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles';
%cycleDataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/';
showPlots=1;
%showPlots=0;

overwriteData=0;

useAdjustedCycleTimes=1;

trackLength=2.5;

useOmarAsymAdjustment=1;

if(useOmarAsymAdjustment)
    adjustedCycleDataDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData/lfpOmarAsymAdjustedCycleInfos';

else
    adjustedCycleDataDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData/lfpAdjustedCycleInfos';

end

    
bySessionSpikeDataDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData/bySessionSpikeData';
touchDir(bySessionSpikeDataDir)

fieldsPerSessionImagesDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3Images/fieldsPerSession';
touchDir(fieldsPerSessionImagesDir)

perFieldDataDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData/perFieldPerCycleData_Dec16_2021';
touchDir(perFieldDataDir)

processedDataDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData';
unitDataDir=fullfile(processedDataDir,'unitSpikeInfoStructs');
dataDir=unitDataDir;

if(~exist('spikeDataPerField','var'))
    spikeDataPerField=load('allUnitInFieldGoodSpeedSpikePhasesAndTimes.mat');
end

spikeDataUnitPathsPerField=spikeDataPerField.unitInfoPathPerField;

sqlDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/sqlQueries';

%Cicero_09172014 has only 1 strongly precessing cell used
%sessionNames={'Achilles_10252013','Achilles_11012013','Buddy_06272013', 'Cicero_09102014','Cicero_09172014','Gatsby_08282013'};
%sessionNames={'Achilles_10252013','Achilles_11012013','Buddy_06272013', 'Cicero_09102014','Gatsby_08282013'};

sqlSessionFiles=getRegexFileNames(sqlDir,'ec*csv');

numSessions=length(sqlSessionFiles);
sessionNames=cell(numSessions,1);

for si=1:numSessions
    sessionNames{si}=getSubstrBtwnSubstrs(sqlSessionFiles{si},'','_sql');
end

halfLineHeight=0.01; %raster line height
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%loop through each session and all corresponding unit-fields
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

startSI=1;
    setTightSubplots
    
    %startSI=37;
    
    
    
    startSI=1;
    
    startSI=78; %ec016.233
    
for si=startSI:length(sessionNames)
    currSessName=sessionNames{si};
    clearvars -except overwriteData trackLength useOmarAsymAdjustment useAdjustedCycleTimes adjustedCycleDataDir spikeDataPerField spikeDataUnitPathsPerField bySessionSpikeDataDir fieldsPerSessionImagesDir perFieldDataDir si currSessName sessionNames showPlots dataDir cycleDataDir
    setTightSubplots
    totalFieldCount=spikeDataPerField.totalFieldCount;

        %{
    currSessMasterUnitData={};
    for fi=1:length(currSessUnitDataPaths)
        currData=load(currSessUnitDataPaths{fi});
        currSessMasterUnitData{fi}=currData;
    end
        %}
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %collect data on units that fire in the same session and direction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sessionAllFieldDataFilePath=fullfile(bySessionSpikeDataDir,sprintf('NEW%sSpikeDataByUnitField.mat',currSessName));
    if(false && exist(sessionAllFieldDataFilePath,'file'))
        load(sessionAllFieldDataFilePath)
    else
        
    %[minTime, maxTime,numLapsPerDir]=getSessionTimeBounds(sessionNames{si});
    %maxTime=minTime+10;
    
    %[leftwardRefUnitData,rightwardRefUnitData]=getRefUnitData(currSessName);
    
    currSessUnitDataPaths=getRegexFilePaths(dataDir,sprintf('%s*Unit*.mat',currSessName));
    
    
        numUnitsInSess=length(currSessUnitDataPaths);
        
        numUnitFieldsInSess=0;
        
        numRightFieldsPerUnitInSession=zeros(numUnitsInSess,1);
        numLeftFieldsPerUnitInSession=zeros(numUnitsInSess,1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %UNIQUE FIELD ID CODE: unitIDinSession*1000+di*100+directionalfieldNum
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        unitPathsPerFieldID={};
        allFieldIDsThisSession=[];
        allSpikeDataFieldIdxesThisSession=[];
        for ui=1:numUnitsInSess
            currUnitFilePath=currSessUnitDataPaths{ui};
            
            spikeDataCandidateIdxes=zeros(totalFieldCount,1);
            for j=1:length(spikeDataUnitPathsPerField)
                
                if(strcmp(currUnitFilePath,spikeDataUnitPathsPerField{j}))
                    
                    spikeDataCandidateIdxes(j)=1;

                end
            end
        
            spikeDataCandidateFieldIdxList=find(spikeDataCandidateIdxes);
            
            fieldIDsInUnit=[];
            spikeDataFieldIdxes=[];
            currUnitData=load(currUnitFilePath);
            unitIDinSession=ui;
            
            numRightFields=length(currUnitData.rightwardFieldStartEndM)/2;
            numLeftFields=length(currUnitData.leftwardFieldStartEndM)/2;
            
            numRightFieldsPerUnitInSession(unitIDinSession)=numRightFields;
            numLeftFieldsPerUnitInSession(unitIDinSession)=numLeftFields;
            
            fieldCountInUnit=0;
            for di=1:2
                if(di==1)
                    currDirNumFields=numRightFields;
                elseif(di==2)
                    currDirNumFields=numLeftFields;
                end
                
                for dfi=1:currDirNumFields
                    fieldIDperUnitPerDir=1000*unitIDinSession+100*di+dfi;
                    fieldIDsInUnit=[fieldIDsInUnit; fieldIDperUnitPerDir];
                    fieldCountInUnit=fieldCountInUnit+1;
                    
                    
                   spikeDataFieldIdxes=[spikeDataFieldIdxes; spikeDataCandidateFieldIdxList(fieldCountInUnit)];
                    unitPathsPerFieldID{fieldIDperUnitPerDir}=currUnitFilePath;
                end
            end

            if(overwriteData)
                save(currUnitFilePath,'unitIDinSession','fieldIDsInUnit','spikeDataFieldIdxes','-append')
            end
             allFieldIDsThisSession=[allFieldIDsThisSession; fieldIDsInUnit(:)];
            allSpikeDataFieldIdxesThisSession=[allSpikeDataFieldIdxesThisSession;spikeDataFieldIdxes(:)];
        end
        
        totalNumRightFieldsInSession=sum(numRightFieldsPerUnitInSession(:));
        totalNumLeftFieldsInSession=sum(numLeftFieldsPerUnitInSession(:));
        
        
        totalNumFieldsInSession=length(allFieldIDsThisSession);
        
        if(totalNumFieldsInSession==0)
            continue
        end
        
        %maxNumFieldsForStorage=max([totalNumRightFieldsInSession totalNumLeftFieldsInSession]);
        
            currSessMasterAllFieldSpikeData=cell(totalNumFieldsInSession,1);
            currSessMasterFieldMiddleData=cell(totalNumFieldsInSession,1);
            
           currSessMasterSpikeNonAdjustedPhases=cell(totalNumFieldsInSession,1);
           currSessMasterSpikePositions=cell(totalNumFieldsInSession,1);
            currSessMasterSpikePositionsFieldFrac=cell(totalNumFieldsInSession,1);
                            
         
                            currSessMasterEnterFieldTimeData=cell(totalNumFieldsInSession,1);
                            currSessMasterExitFieldTimeData=cell(totalNumFieldsInSession,1);
                 
                             
                   
            %currSessMasterEnterFieldTimeData={};
            %currSessMasterFirstSpikeData={};

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%for each direction, loop through all units and collect positional info
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %currSessMasterAllFieldSpikeData=cell(maxNumFieldsForStorage,2);
        %currSessMasterFieldMiddleData=cell(maxNumFieldsForStorage,2);
        %currSessMasterEnterFieldTimeData=cell(maxNumFieldsForStorage,2);
        %currSessMasterExitFieldTimeData=cell(maxNumFieldsForStorage,2);

        fieldCenterPerField=NaN(totalNumFieldsInSession,1);
        fieldWidthPerField=NaN(totalNumFieldsInSession,1);
        fieldStartPerField=NaN(totalNumFieldsInSession,1);
        fieldEndPerField=NaN(totalNumFieldsInSession,1);
        
                for fi=1:totalNumFieldsInSession
                            currFieldID=allFieldIDsThisSession(fi);
                            
                            currFieldSpikeDataIdx=allSpikeDataFieldIdxesThisSession(fi);
                            [ui,di,dfi]=recoverFieldInfoFromID(currFieldID);
                            
                            %currUnitData=currSessMasterUnitData{fi};
                            currUnitData=load(unitPathsPerFieldID{currFieldID});
                            
                       
                             if(di==1)
                                currDir='right';
                                currDirStr='rightward';
                                 thisUnitSpikePhases=currUnitData.unitInfo.rightSpikePhases;
                                 thisUnitSpikePositions=currUnitData.unitInfo.rightSpikePositions;
                            else
                                currDir='left';
                                currDirStr='leftward';
                                thisUnitSpikePhases=currUnitData.unitInfo.leftSpikePhases;
                                thisUnitSpikePositions=currUnitData.unitInfo.leftSpikePositions;
                             end
                            
                            allTimesStr=(sprintf('%sSpikeTimes',currDir));
                            %firstTimesStr=(sprintf('%sSpikeTimesOnlyRefCycleFirsts',currDir));
                            
                            if(~isfield(currUnitData,'directionSpecificStats'))
                                continue
                            end
                            if(isfield(currUnitData.directionSpecificStats,currDirStr))
                            %fieldMiddlePosPerTime=(currUnitData.directionSpecificStats.(currDirStr).fieldExitPosPerTime + ...
                             %   currUnitData.directionSpecificStats.(currDirStr).fieldEntryPosPerTime)/2;
                                fieldNumPerTime=currUnitData.directionSpecificStats.(currDirStr).fieldNumPerTime;
                                fieldEntryTimePerTime=currUnitData.directionSpecificStats.(currDirStr).fieldEntryTimePerTime;
                                fieldExitTimePerTime=currUnitData.directionSpecificStats.(currDirStr).fieldExitTimePerTime;
                                
                                numFields=max(fieldNumPerTime);
                                %fieldCenters=(currUnitData.(sprintf('manualFieldStartsM%s',currDirStr))+currUnitData.((sprintf('manualFieldEndsM%s',currDirStr))))/2;
                                fieldStartStops=currUnitData.(sprintf('%sFieldStartEndM',currDirStr));
                                
                                thisFieldPosStart=fieldStartStops(2*(dfi-1)+1);
                                thisFieldPosEnd=fieldStartStops(2*(dfi-1)+2);
                                
                                thisFieldCenter=(thisFieldPosStart+thisFieldPosEnd)/2;
                                thisFieldWidth=abs(thisFieldPosStart-thisFieldPosEnd);
                                
                                if(di==1)
                                    inFieldSpikeIdxes=thisUnitSpikePositions>=thisFieldPosStart & thisUnitSpikePositions <=thisFieldPosEnd;
                                elseif(di==2)
                                   inFieldSpikeIdxes=thisUnitSpikePositions<=thisFieldPosStart & thisUnitSpikePositions >=thisFieldPosEnd;
                                end
                                
                                inFieldSpikePositions=thisUnitSpikePositions(inFieldSpikeIdxes);
                                
                                
                                inFieldSpikePositionsFieldFrac=abs(thisUnitSpikePositions(inFieldSpikeIdxes)-thisFieldPosStart)/thisFieldWidth;
                                
                                inFieldSpikePhases=thisUnitSpikePhases(inFieldSpikeIdxes);
                                
                                inFieldSpikePhasesInSpikeData=spikeDataPerField.spikePhasesPerField{currFieldSpikeDataIdx};
                                
                                spikeDataIdxCorrect=1;
                                numDiff=0;
                                for sdi=1:length(inFieldSpikePhasesInSpikeData)
                                    if(min(abs(inFieldSpikePhasesInSpikeData(sdi)-inFieldSpikePhases))>0)
                                        numDiff=numDiff+1;
                                        %spikeDataIdxCorrect=0;
                                    end
                                end
                                
                                if(numDiff<length(inFieldSpikePhasesInSpikeData)/3)
                                    spikeDataIdxCorrect=1;
                                end
                                if(spikeDataIdxCorrect==0)
                                    disp('')
                                end
                                
                                %assert(spikeDataIdxCorrect==1)
                                
                                numBounds=length(fieldStartStops);
                                fieldCenters=NaN(numBounds/2,1);
                                for bi=1:2:(numBounds-1)
                                    currFieldIdx=(bi+1)/2;
                                    fieldCenters(currFieldIdx)=(fieldStartStops(bi)+fieldStartStops(bi+1))/2;
                                end
                                
                                 fieldMiddlePosPerTime=NaN(size(fieldNumPerTime));
                                for ti=1:length(fieldNumPerTime)
                                    if(~isnan(fieldNumPerTime(ti)))
                                     fieldMiddlePosPerTime(ti)=fieldCenters(fieldNumPerTime(ti));
                                    end
                                end
                            else
                                currSessMasterAllFieldSpikeData{fi}=NaN;
                                currSessMasterFieldMiddleData{fi}=NaN;
                                currSessMasterEnterFieldTimeData{fi}=NaN;
                                currSessMasterExitFieldTimeData{fi}=NaN;
                                continue
                            end
                            
                            if(isfield(currUnitData,allTimesStr))
                                currSessMasterAllFieldSpikeData{fi}=currUnitData.(allTimesStr)(inFieldSpikeIdxes); %ONLY IN FIELD SPIKE TIMES
                            else
                                currSessMasterAllFieldSpikeData{fi}=NaN;
                            end
                            
                            currSessMasterSpikeNonAdjustedPhases{fi}=inFieldSpikePhases;
                            currSessMasterSpikePositions{fi}=inFieldSpikePositions;
                            currSessMasterSpikePositionsFieldFrac{fi}=inFieldSpikePositionsFieldFrac;
                            
                            currSessMasterFieldMiddleData{fi}=fieldMiddlePosPerTime;
                            currSessMasterEnterFieldTimeData{fi}=fieldEntryTimePerTime;
                            currSessMasterExitFieldTimeData{fi}=fieldExitTimePerTime;
                            
                            
                            fieldCenterPerField(fi)=thisFieldCenter;
                            fieldWidthPerField(fi)=thisFieldWidth;
                            fieldStartPerField(fi)=thisFieldPosStart;
                            fieldEndPerField(fi)=thisFieldPosEnd;
                            
                            %spikeDataFieldIdxPerField(fi)=
                            %{
                            if(isfield(currUnitData,firstTimesStr))
                                currSessMasterFirstSpikeData{fi,di}=currUnitData.(firstTimesStr);
                            else
                                currSessMasterFirstSpikeData{fi,di}=NaN;
                            end
                            %}
                            %figure; imshow(currUnitData.unitInfo.comboImageFilePath)
                end
                %toc
   
            %save(fullfile(bySessionSpikeDataDir,sprintf('%sSpikeData.mat',currSessName)))
            if(overwriteData)
                save(sessionAllFieldDataFilePath)
            end
        end


        %{
refUnitData=leftwardRefUnitData;
    if(~isstruct(refUnitData))
            continue
        end
        %}
        
    avgThetaAmpPerUnit=NaN(numUnitsInSess,1);
    for ui=1:numUnitsInSess

        unitData=load(currSessUnitDataPaths{ui}); 
        currUnitCycleInfo=load(unitData.refChCycleInfoPath); 
        
        if(isfield(currUnitCycleInfo,'ampPerCycle'))
            ampPerCycle=currUnitCycleInfo.ampPerCycle;
        else
            ampPerCycle=(currUnitCycleInfo.cycleMaxAmp - currUnitCycleInfo.cycleMinAmp); 
            if(overwriteData)
                save(unitData.refChCycleInfoPath,'ampPerCycle','-append')
            end
        end
        avgAmp=nanmean(ampPerCycle);
        avgThetaAmpPerUnit(ui)=avgAmp;
    end
    
    [~,bestRefUnitID]=max(avgThetaAmpPerUnit);
    
    refUnitData=load(currSessUnitDataPaths{bestRefUnitID}); 
    
    

    if(useAdjustedCycleTimes)
            adjustedCycleInfoFilePath=getRegexFilePath(adjustedCycleDataDir,sprintf('*_%s_adjCycleLFPstruct.mat', currSessName));
            cycleInfoStruct=load(adjustedCycleInfoFilePath);
          refChStr=lower(getSubstrBtwnSubstrs(refUnitData.refChCycleInfoPath,sprintf('%s_',currSessName),'_6-12_minmax.mat'));

          refChCycleInfo=cycleInfoStruct.adjustedCycleInfoPerPyrCh.(refChStr);
    else
        refChCycleInfo=load(refUnitData.refChCycleInfoPath);
    end

    %cycleIsDuringMaze=refChCycleInfo.cycleMinTimes>=minTime & refChCycleInfo.cycleMinTimes<maxTime;

    if(useAdjustedCycleTimes)
        if(useOmarAsymAdjustment)
            cycleInMazeStartTimes=refChCycleInfo.fastStartTimes;
            cycleInMazeEndTimes=refChCycleInfo.fastEndTimes;
        else
            cycleInMazeStartTimes=refChCycleInfo.closestFastMaxTimes(1:(end-1));
            cycleInMazeEndTimes=refChCycleInfo.closestFastMaxTimes(2:end);
        end
    else
        cycleInMazeStartTimes=refChCycleInfo.cycleMaxTimes(1:(end-1));
        cycleInMazeEndTimes=refChCycleInfo.cycleMaxTimes(2:end);
    end

    %cycleInMazeStartTimes=cycleStartTimes(cycleIsDuringMaze);
    %cycleInMazeEndTimes=cycleEndTimes(cycleIsDuringMaze);

    numCycles=length(cycleInMazeStartTimes);
    
    posInfo=load(refUnitData.unitInfo.positionInfoForSessionSavePath);

    disp('')
    
    %filledSignedSpeedPerTime=refUnitData.unitInfo.filledSignedSpeedPerTime;
    filledSignedSpeedPerTime=refUnitData.unitInfo.speedPerTimeStep;
    positionTimeAxis=posInfo.positionTimeAxisSec;
    
    minSpeed=0.05;
    
     allSpikeTimesByCycleAndFieldID={};
     allSpikeTimesByCyclePerField=cell(numCycles,totalNumFieldsInSession);
     
     peakSpikePhasePerCyclePerField=NaN(numCycles,totalNumFieldsInSession);
     avgSpikeTimesPerCyclePerField=NaN(numCycles,totalNumFieldsInSession);
     avgSpikeFieldFracPosPerCyclePerField=NaN(numCycles,totalNumFieldsInSession);
     
     fieldCenterPerCyclePerField=NaN(numCycles,totalNumFieldsInSession);
     fieldEntryTimePerCyclePerField=NaN(numCycles,totalNumFieldsInSession);
     fieldExitTimePerCyclePerField=NaN(numCycles,totalNumFieldsInSession);
     speedPerCycle=NaN(numCycles,1);
     timeOfCycleTrough=NaN(numCycles,1);
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %loop through each theta cycle 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic
 

    for ci=1:numCycles
        
        currCycleStartTime=cycleInMazeStartTimes(ci);
        currCycleEndTime=cycleInMazeEndTimes(ci);
        cycleCenterTime=(currCycleStartTime+currCycleEndTime)/2;
        
        currCycleDuration=currCycleEndTime-currCycleStartTime;
        
        currSpeed=interp1(positionTimeAxis,filledSignedSpeedPerTime,cycleCenterTime);
        %currPos=interp1(positionTimeAxis, posAlongTrackPerTimeM,cycleCenterTime);
        
        speedPerCycle(ci)=abs(currSpeed);
        timeOfCycleTrough(ci)=cycleCenterTime;
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %SPEED THRESHOLD
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if(currSpeed>minSpeed)
            di=1; %rightward movement
        elseif(currSpeed<-minSpeed)
            di=2;%leftward movement
        else
            %allSpikeTimesByCycleAndFieldID{ci,fi}=NaN;
            continue
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %find direction during this cycle
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %for di=1:2
        if(di==1)
            currDirStr='rightward';
            currDir='right';
        else
            currDirStr='leftward';
            currDir='left';
        end
        
        if(isempty(currSessMasterAllFieldSpikeData))
            continue
        end

    %positionTimeAxis=posInfo.positionTimeAxis;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %loop through each field in this session, gather spikes in
        %current cycle
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for fi=1:totalNumFieldsInSession
            currFieldAllSpikeTimes=currSessMasterAllFieldSpikeData{fi};
            currFieldAllFieldMiddlePerTime=currSessMasterFieldMiddleData{fi};
            
            currFieldEntryTimePerTime=currSessMasterEnterFieldTimeData{fi};
            currFieldExitTimePerTime=currSessMasterExitFieldTimeData{fi};
            
           currFieldSpikeFieldFracPositions=currSessMasterSpikePositionsFieldFrac{fi};

            
            currFieldID=allFieldIDsThisSession(fi);
             [ui,di,dfi]=recoverFieldInfoFromID(currFieldID);
             
            if(isempty(currFieldAllFieldMiddlePerTime))
                continue
            end
            if(isnan(currFieldAllFieldMiddlePerTime))
                continue
            end
            if(isempty(currFieldAllSpikeTimes))
                continue
            end
            if(isnan(currFieldAllSpikeTimes))
                continue
            end
            
            inCycleTimeIdxes=positionTimeAxis>=currCycleStartTime & positionTimeAxis<currCycleEndTime;
            currCycleFieldMid=currFieldAllFieldMiddlePerTime(inCycleTimeIdxes);
            
            currFieldCycleFieldEntryTimes=currFieldEntryTimePerTime(inCycleTimeIdxes);
            currFieldCycleFieldExitTimes=currFieldExitTimePerTime(inCycleTimeIdxes);
            
            inCycleSpikeIdxes=currFieldAllSpikeTimes>=currCycleStartTime & currFieldAllSpikeTimes<currCycleEndTime;
            currFieldAndCycleSpikeTimes=currFieldAllSpikeTimes(inCycleSpikeIdxes);
            
            
            currFieldAndCycleSpikeFieldFracPositions=currFieldSpikeFieldFracPositions(inCycleSpikeIdxes);
            
            currUnitDirCycleSpikePhases=360*(currFieldAndCycleSpikeTimes-currCycleStartTime)/currCycleDuration;
            
             if(~isempty(currUnitDirCycleSpikePhases))
                kernelDistrPhaseObj = fitdist(currUnitDirCycleSpikePhases(:),'Kernel','BandWidth',30); %gaussian kernel estimation, 30 degree
                %phaseAxis=linspace(0,360,361);
                phaseAxis=linspace(0,360,1000);
                kernelDistrPhase=pdf(kernelDistrPhaseObj,phaseAxis);

                maxProb=max(kernelDistrPhase);
                maxIdxes=find(kernelDistrPhase==maxProb);

                peakSpikePhasePerCyclePerField(ci,fi)=nanmean(phaseAxis(maxIdxes));
            else
                peakSpikePhasePerCyclePerField(ci,fi)=NaN;
            end

            allSpikeTimesByCyclePerField{ci,fi}=currFieldAndCycleSpikeTimes-currCycleStartTime;
            avgSpikeTimesPerCyclePerField(ci,fi)=nanmean(currFieldAndCycleSpikeTimes-currCycleStartTime);
            avgSpikeFieldFracPosPerCyclePerField(ci,fi)=nanmean(currFieldAndCycleSpikeFieldFracPositions);
            
             %avgSpikeTimesByCycleAndUnit(ci,fi,di)=nanmean(min(currUnitDirCycleSpikeTimes))-currCycleStartTime;

            fieldCenterPerCyclePerField(ci,fi)=nanmean(currCycleFieldMid);
            
            fieldEntryTimePerCyclePerField(ci,fi)=nanmean(currFieldCycleFieldEntryTimes);
            fieldExitTimePerCyclePerField(ci,fi)=nanmean(currFieldCycleFieldExitTimes);
        end
        
        if(mod(ci,1000)==0)
            disp(ci/numCycles)
        end
    end
    
   
    toc
    
    minNumCellsPerSeq=5;
    minNumCellsPerSeq=1;
    onlyInFieldSpiking=1;
    
    if(onlyInFieldSpiking)
        avgSpikeTimesPerCyclePerField(isnan(fieldCenterPerCyclePerField))=NaN;
        avgSpikeFieldFracPosPerCyclePerField(isnan(fieldCenterPerCyclePerField))=NaN;
        peakSpikePhasePerCyclePerField(isnan(fieldCenterPerCyclePerField))=NaN;
    end
    
    
    %plot raster across session per lap...?
    if(showPlots)
        J=customcolormap_preset('pasteljet',250);
        
        %fieldColors=jet(250); %250 cm
        fieldColors=J; %250 cm
        
        
        %fR=figure(11);
        %fL=figure(12);
        
        figure
s1=subplot(2,1,1)
s2=subplot(2,1,2)
        for fi=1:totalNumFieldsInSession
             currFieldID=allFieldIDsThisSession(fi);
             [ui,di,dfi]=recoverFieldInfoFromID(currFieldID);
             
             if(di==1)
                currFieldCenterPos=round(fieldCenterPerField(fi)*100);
                currFieldSpikePositions=currSessMasterSpikePositions{fi};
             else
                 currFieldCenterPos=round((trackLength-fieldCenterPerField(fi))*100);
                 currFieldSpikePositions=trackLength-currSessMasterSpikePositions{fi};
             end
             
             
             currFieldSpikePhases=currSessMasterSpikeNonAdjustedPhases{fi};
             
             if(di==1)
                 %figure(fR)
                 
                 axes(s1)
                 hold on
                 dirStr='rightward';
             else
                 %figure(fL)
                 axes(s2)
                 hold on
                 dirStr='leftward';
             end
             
            % plot(currFieldSpikePositions,currFieldSpikePhases,'.','MarkerSize',20,'Color',fieldColors(currFieldCenterPos,:))
             plot(currFieldSpikePositions,currFieldSpikePhases,'.','MarkerSize',12,'Color',fieldColors(currFieldCenterPos,:))
             drawnow
             box off

             hold on
             
             title(sprintf('%s, all %s fields in session',currSessName,dirStr))
             ylabel('Theta phase (deg)')
             ylim([0 360])
             
             xlabel('Position (m)')
             xlim([0 2.5])
             colormap(jet)
              cb=colorbar;
              caxis([0 2.5])
             ylabel(cb,'Field center (m)')
            
             
        end

        maxFig
        setFigFontTo(24)
       

     
    end
     allFieldsInSessionImgPath=fullfile(fieldsPerSessionImagesDir,sprintf('%s_AllFields.png',currSessName));

     
     if(showPlots)
        saveas(gcf,allFieldsInSessionImgPath)
        saveEPS(fullfile(fieldsPerSessionImagesDir,sprintf('%s_AllFields.eps',currSessName)))
     end
   
    %save(sprintf('perUnitPerCycleTimingData%s.mat',currSessName),'allFieldIDsThisSession','fieldCenterPerCyclePerField','fieldEntryTimePerCyclePerField','fieldExitTimePerCyclePerField','avgSpikeTimesPerCyclePerField','numCycles','onlyInFieldSpiking','speedPerCycle','sortedMeanCenterIdxPerDi','validCycleIdxesPerDi','timeOfCycleTrough','peakSpikePhasePerCyclePerField')
      if(overwriteData)
    save(fullfile(perFieldDataDir,sprintf('perFieldPerCycleTimingData%s.mat',currSessName)),...
           'allFieldsInSessionImgPath','fieldStartPerField','fieldEndPerField','allFieldIDsThisSession','allSpikeDataFieldIdxesThisSession',...
           'totalNumFieldsInSession','fieldCenterPerField','fieldWidthPerField','currSessMasterSpikePositions', ...
           'currSessMasterSpikePositionsFieldFrac','currSessMasterSpikeNonAdjustedPhases','fieldCenterPerCyclePerField',...
           'fieldEntryTimePerCyclePerField','fieldExitTimePerCyclePerField','avgSpikeTimesPerCyclePerField',...
           'avgSpikeFieldFracPosPerCyclePerField','numCycles','onlyInFieldSpiking','speedPerCycle','timeOfCycleTrough',...
           'peakSpikePhasePerCyclePerField','allSpikeTimesByCyclePerField','useAdjustedCycleTimes')
       %save(sprintf('perFieldPerCycleTimingData%s.mat',currSessName),'fieldStartPerField','fieldEndPerField','-append')
      end
    close all
end

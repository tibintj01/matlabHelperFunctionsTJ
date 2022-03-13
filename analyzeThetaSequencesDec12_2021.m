close all; clear all; clc
%dataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles';

processedDataDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData';
unitDataDir=fullfile(processedDataDir,'unitSpikeInfoStructs');
dataDir=unitDataDir;

saveThetaSeqImgDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3Images/individualThetaSequencesSimilarSizeMinTrajLength';
touchDir(saveThetaSeqImgDir)

posInfoDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData/positionTimeInfos';
%perCycleDataDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData/perCycleUnitSeqData';
perCycleDataDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData/perFieldPerCycleData';

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
minTrajLengthWidthRatio=1/3;

for si=startSi:length(sessionNames) 
    
    currSessName=sessionNames{si};
    currSeqDataFilePath=fullfile(perCycleDataDir,sprintf('perFieldPerCycleTimingData%s.mat',currSessName));
    
    currPosInfoFilePath=getFilePathsRegex(posInfoDir,sprintf('*%s*Position*.mat',currSessName));
    
    representativeUnitFilePaths=getFilePathsRegex(unitDataDir,sprintf('*%s*.mat',currSessName));

    if(isempty(representativeUnitFilePaths))
        continue
    end
    
    loadRepresentativeUnitFilesForSess
    
    
    leftwardLapStartTimes=representativeLeftwardUnitInfo.lapStartTimesPerDir.leftward;
    leftwardLapStopTimes=representativeLeftwardUnitInfo.lapStopTimesPerDir.leftward;
    
    rightwardLapStartTimes=representativeRightwardUnitInfo.lapStartTimesPerDir.rightward;
    rightwardLapStopTimes=representativeRightwardUnitInfo.lapStopTimesPerDir.rightward;
    
    currSessionPosInfo=load(currPosInfoFilePath);
    
    currThetaSeqData=load(currSeqDataFilePath);
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
          
          %approxLapTime=20; %sec
    for ci=1:numCycles
        
        currCycleTime=timeOfCycleTrough(ci);
        
        approxCurrLapStartTime=currCycleTime-approxLapTime/2;
         approxCurrLapEndTime=currCycleTime+approxLapTime/2;
         
         currLapPosTimeIdxes=posTimeAxis>approxCurrLapStartTime & posTimeAxis<approxCurrLapEndTime;
        
        currCyclePeakPhasePerField=peakSpikePhasePerCyclePerField(ci,:);
        
        %at least n fields in sequence
        if(sum(~isnan(currCyclePeakPhasePerField))<minNumFieldsInSeq)
            continue
        end
        
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
            currCycleFieldCenterPerFieldInSpaceOrder=fieldCenterPerField(rightwardSortedFieldCenterIdxes);
            ratCurrCyclePos=interp1(posTimeAxis,posPerTime,currCycleTime);
            
            lapStartTimes=rightwardLapStartTimes;
            lapStopTimes=rightwardLapStopTimes;
            
            currDirSortedFieldCenterIdxes=rightwardSortedFieldCenterIdxes;

        elseif(allFieldsLeftward)
            currCyclePeakPhasePerFieldInSpaceOrder=currCyclePeakPhasePerField(leftwardSortedFieldCenterIdxes);
            currCycleFieldCenterPerFieldInSpaceOrder=fieldCenterPerField(leftwardSortedFieldCenterIdxes);
            
            ratCurrCyclePos=interp1(posTimeAxis,maxTrackLength-posPerTime,currCycleTime);
            
            lapStartTimes=leftwardLapStartTimes;
            lapStopTimes=leftwardLapStopTimes;
            
            currDirSortedFieldCenterIdxes=leftwardSortedFieldCenterIdxes;
        end
        
        fieldStartPerFieldInSpaceOrder=fieldStartPerField(currDirSortedFieldCenterIdxes);
        fieldEndPerFieldInSpaceOrder=fieldEndPerField(currDirSortedFieldCenterIdxes);
        
        currCycleActiveFieldOrderedIdxes=~isnan(currCyclePeakPhasePerFieldInSpaceOrder);
        
        currCyclePhaseSeq=currCyclePeakPhasePerFieldInSpaceOrder(currCycleActiveFieldOrderedIdxes);
        

        if(allFieldsLeftward)
            currCycleFieldCenterSeq=maxTrackLength-currCycleFieldCenterPerFieldInSpaceOrder(currCycleActiveFieldOrderedIdxes); %1-distances for leftward runs
            currCycleFieldStartSeq=maxTrackLength-fieldStartPerFieldInSpaceOrder(currCycleActiveFieldOrderedIdxes);
            currCycleFieldEndSeq=maxTrackLength-fieldEndPerFieldInSpaceOrder(currCycleActiveFieldOrderedIdxes);
           
        else
             currCycleFieldCenterSeq=currCycleFieldCenterPerFieldInSpaceOrder(currCycleActiveFieldOrderedIdxes);
            currCycleFieldStartSeq=fieldStartPerFieldInSpaceOrder(currCycleActiveFieldOrderedIdxes);
            currCycleFieldEndSeq=fieldEndPerFieldInSpaceOrder(currCycleActiveFieldOrderedIdxes);
        end
        
        currCycleFieldWidthSeq=currCycleFieldEndSeq-currCycleFieldStartSeq;
        
        %currCycleFieldCenterSeqAroundCenter=currCycleFieldCenterSeq-nanmean(currCycleFieldCenterSeq);
        
        angDiffPerField=angdiffDeg(currCyclePhaseSeq);
        
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
        
        if(enforceMinTrajectoryLength)
            if(trajLengthFieldWidthRatio<minTrajLengthWidthRatio)
                continue
            end
            
        end
        
        if(showPlots)
             fH=figure;
             
            %plot(currCyclePhaseSeq,ones(size(currCyclePhaseSeq)),'k.')
            currCycleFieldCount=length(currCyclePhaseSeq);
            
            subplot(1,2,1)
           distColors=jet(length(currCycleFieldCenterSeq))
           for ai=1:length(currCycleFieldCenterSeq)
               %plot(currCycleFieldCenterSeq(ai),currCyclePhaseSeq(ai),'.','Color',distColors(ai,:),'MarkerSize',50)
               plot(currCyclePhaseSeq(ai),currCycleFieldCenterSeq(ai),'.','Color',distColors(ai,:),'MarkerSize',50)
               hold on
               %plot([currCycleFieldStartSeq(ai) currCycleFieldEndSeq(ai)],[currCyclePhaseSeq(ai) currCyclePhaseSeq(ai)],'-','Color',distColors(ai,:),'LineWidth',5)
               plot([currCyclePhaseSeq(ai) currCyclePhaseSeq(ai)],[currCycleFieldStartSeq(ai) currCycleFieldEndSeq(ai)],'-','Color',distColors(ai,:),'LineWidth',5)
 
               %vline(
           end
           colormap(jet)
           cb=colorbar
            ylabel(cb,'Field order')
            caxis([1 length(currCycleFieldCenterSeq)])

            %hold on
            %xlim([0 maxTrackLength])
            distRange=range([min(currCycleFieldStartSeq) max(currCycleFieldEndSeq)]);
            
            %{
            xlim([min(currCycleFieldStartSeq)-distRange*0.2 max(currCycleFieldEndSeq)+distRange*0.2])
            xlabel('Distance along track (m)')
            %ylim([-0.5 0.5])
            ylim([0 360])
            ylabel('Theta phase (deg)')
            %}
            ylim([min(currCycleFieldStartSeq)-distRange*0.2 max(currCycleFieldEndSeq)+distRange*0.2])
            ylabel('Distance along track (m)')
            %ylim([-0.5 0.5])
            xlim([0 360])
            xlabel('Theta phase (deg)')
            
             
            %pL=vline(ratCurrCyclePos,'k--',3)
            pL=hline(ratCurrCyclePos,'k--',3)
            hold off
            
            title(sprintf('%s, single theta sequence, cycle %d',currSessName, ci))
            %legend([pL pF],{'current position','active field'},'Location','northwest')%'eastoutside')
            legend([pL],{'current position'},'Location','northwest')%'eastoutside')

            legend boxoff
            axis square
            box off
            
            subplot(1,2,2)
            %plot(posPerTime(currLapPosTimeIdxes),speedPerTime(currLapPosTimeIdxes),'b-')
            
            timeVals=posTimeAxis(currLapPosTimeIdxes)-currCycleTime;
            if(allFieldsLeftward)
                %plotWithColorVec(maxTrackLength - posPerTime(currLapPosTimeIdxes),speedPerTime(currLapPosTimeIdxes),timeVals,'jet',5)
                plotWithColorVec(speedPerTime(currLapPosTimeIdxes),maxTrackLength - posPerTime(currLapPosTimeIdxes),timeVals,'jet',5)
            else
               %plotWithColorVec(posPerTime(currLapPosTimeIdxes),speedPerTime(currLapPosTimeIdxes),timeVals,'jet',5)
                plotWithColorVec(speedPerTime(currLapPosTimeIdxes),posPerTime(currLapPosTimeIdxes),timeVals,'jet',5)

            end
            
            cb=colorbar
            ylabel(cb,'Time from current cycle (sec)')
            
            %{
            xlim([min(currCycleFieldStartSeq)-distRange*0.2 max(currCycleFieldEndSeq)+distRange*0.2])
            xlabel('Distance along track (m)')
              ylim([0 1])
            ylabel('Running speed (m/s)')
            %}
            
             ylim([min(currCycleFieldStartSeq)-distRange*0.2 max(currCycleFieldEndSeq)+distRange*0.2])
            ylabel('Distance along track (m)')
              xlim([0 1])
            xlabel('Running speed (m/s)')
            hold on
            
            %pL=vline(ratCurrCyclePos,'k--',3)
            pL=hline(ratCurrCyclePos,'k--',3)
            
            legend(pL,'current position','Location','northwest')
            legend boxoff
             title(sprintf('running speed over time'))
             
            axis square
            hold off
            
            setFigFontTo(18)
            %maxFigHalfHalfWidth
            maxFig
            box off
            saveas(gcf,fullfile(saveThetaSeqImgDir,sprintf('%s_ThetaCycle%d_%dFields.png',currSessName, ci,currCycleFieldCount)))
            %drawnow
            %pause(0.3)
            close all
        end
  
    end
               
end
        
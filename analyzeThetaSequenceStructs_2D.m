close all;  clc
%dataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles';

clearvars -except spikeDataPerField

offsetAndDiffVsSpeedDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3Images/offsetAndDiffVsSpeedIndThetaSeqMin2FieldsCntrlDist';

%exampleThetaSeqVsSpeedDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3Images/exampleThetaSeqVsSpeed_2DthetaSequences_TopVsBottom20PrctileSpeed_TimeInFieldOneDenom_MinTrajLength1Eighth_CntrlForPosition';
%exampleThetaSeqVsSpeedDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3Images/exampleThetaSeqVsSpeed_2DthetaSequences_TopVsBottom25PrctileSpeed_NoVarRestr_TimeInFieldOneDenom_MinTrajLengthStartAndEnd1Eighth_CntrlForPositionTolSixth_WithSpeedSignatures_MaxSpeed1.5_OmarAsymAdjustedThetaWaveform_WEdges_RealTimeToFieldEnd';
%exampleThetaSeqVsSpeedDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3Images/exampleThetaSeqVsSpeed_2DthetaSequences_TopVsBottom25PrctileSpeed_NoVarRestr_TimeInFieldOneDenom_MinTrajLengthStartAndEnd1Eighth_NoCntrlForPosition_WithSpeedSignatures_MaxSpeed1.5_OmarAsymAdjustedThetaWaveform_WEdges_RealTimeToFieldEnd';
exampleThetaSeqVsSpeedDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3Images/exampleThetaSeqVsSpeed_2DthetaSequences_TopVsBottom25PrctileSpeed_NoVarRestr_TimeInFieldOneDenom_NoMinTrajLength_CntrlForPositionTolSixth_WithSpeedSignatures_MaxSpeed2_OmarAsymAdjustedTheta_WEdges_RealTimeToFieldEnd_WidthSeventyFivePrcnt_CompNew';

touchDir(exampleThetaSeqVsSpeedDir)

%nonCircStats=1;
nonCircStats=0;

usePeakPhaseAsSummary=0;

controlForDistanceAcrossLaps=1;
%controlForDistanceAcrossLaps=0;

touchDir(offsetAndDiffVsSpeedDir)
tic
disp('loading spike data per field')
if(~exist('spikeDataPerField','var'))
    spikeDataPerField=load('allUnitInFieldGoodSpeedSpikePhasesAndTimes.mat');
end
toc

maxRunSpeed=1;

minSpeed=0.1;
minSpeed=0.125;
%maxSpeed=0.7;
maxSpeed=0.7;

speedInFieldsPerSec=1;

%noEdges=1;
noEdges=0;

if(speedInFieldsPerSec)
    %minSpeed=0.2;
    %minSpeed=0.3;
    %minSpeed=0.25;
    minSpeed=0.2;
    %maxSpeed=1.2;
     %maxSpeed=1.25;
     %maxSpeed=1.2;
     maxSpeed=1.4;
     %maxSpeed=1.5;
     %maxSpeed=2;
end

%zScoreSpeed=1;
zScoreSpeed=0;

   minZspeed=-2.5;
    maxZspeed=2;
    
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


%setTightSubplots
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
useCircMeasureSpread=1;
useCircMeasureSpread=0;

if(zScoreSpeed)
    minSpeedStdDevAcrossLaps=0.05;
else
    minSpeedStdDevAcrossLaps=0;
end
%maxSpeedStdDevWithinTrav=0.1;
maxSpeedVarWithinTrav=1;
%maxSpeedStdDevWithinTrav=0.1;

maxSpeedVarWithinTrav=0.03; %sem within 3cm/s
maxSpeedVarWithinTrav=0.025; %sem within 2cm/s

%maxSpeedVarWithinTrav=Inf; %sem within 2cm/s


maxSpeedVarWithinTrav=0.125; %std dev
%maxSpeedVarWithinTrav=0.05; %std dev
%maxSpeedVarWithinTrav=0.15; %std dev

%maxSpeedVarWithinTrav=0.2; %std dev
maxSpeedVarWithinTrav=Inf; 

minNumPtsInSpeedCategory=2;
minNumPtsInSpeedCategory=1;

%minSpeedStdDevAcrossLaps=0.15;
minNumOccurences=5;
minNumFieldsInSeq=2;
%minNumFieldsInSeq=4;
showPlots=1;
showPlots=0;
showThetaWaveforms=0;
 showSpeedSortedRaster=1;

%setTightSubplots_Spacious

if(showPlots)
    fH=figure;
end
masterTravSpeed=[];
masterPhaseDiff=[];
masterMeanPhase=[];
masterPhaseDiffSpeedRhos=[];
masterMeanPhaseSpeedRho=[];
masterPhaseDiffSpeedSlope=[];
masterMeanPhaseSpeedSlope=[];

masterPhaseDiffHighMinusLowSpeed=[];
masterMeanPhaseHighMinusLowSpeed=[];
masterMeanDistFracHighMinusLowSpeed=[];
indFieldHighSpeedPrctileCutoff=80;
indFieldLowSpeedPrctileCutoff=20;

indFieldHighSpeedPrctileCutoff=(2/3)*100;
indFieldLowSpeedPrctileCutoff=(1/3)*100;

indFieldHighSpeedPrctileCutoff=(3/4)*100;
indFieldLowSpeedPrctileCutoff=(1/4)*100;

%indFieldHighSpeedPrctileCutoff=90;
%indFieldLowSpeedPrctileCutoff=10;

maxAllowableSlopeMagnitude=720;
maxAllowableSlopeMagnitude=Inf;


allMeanPhaseSpeedRho=[];
allMeanPhaseSpeedPval=[];
allPhaseSepSpeedRho=[];
allPhaseSepSpeedPval=[];

allSessionRunSpeed=[];
allSessionMeanDistFracs=[];
allSessionAvgDistSepFracs=[];
allSessionAvgTimeSepFracs=[];
 allSessionPhaseSep=[];
 allSessionMeanPhase=[];
 allSessionPhaseComp=[];
 
runSpeedsPerFieldSequence=[];
 phaseSepPerFieldSequence=[];
 meanPhasePerFieldSequence=[];
 phaseCompPerFieldSequence=[];
 meanDistPerFieldSequence=[];
 distSepPerFieldSequence=[];
 timeSepPerFieldSequence=[];
 
 fieldSeqCount=0;

 

 
  
 %exampleSessionName='ec013.762';
 %exampleFieldSeqStr='fieldSeq_861_870_864_868';
 
  exampleSessionName='ec013.761';
 exampleFieldSeqStr='fieldSeq_842_831_841';
 
  %fig 5 example
 %ec013.718_fieldSeq_746_743
 exampleSessionName='ec013.718';
 exampleFieldSeqStr='fieldSeq_746_743';
 
 exampleOnly=1;
  exampleOnly=0;

startSi=1;
%startSi=39;
 
for si=startSi:length(sessionNames) 
    si
    currSessName=sessionNames{si};
    currSeqDataFilePath=fullfile(perCycleDataDir,sprintf('perFieldPerCycleTimingData%s.mat',currSessName));
    
    
    if(exampleOnly && ~strcmp(currSessName,exampleSessionName))
        continue
    end
    
    %missing seq files?
    if(exist(currSeqDataFilePath,'file'))
        currThetaSeqData=load(currSeqDataFilePath);
    else
        continue
    end
    
    %{
    if(~strcmp(currSessName,'ec013.626'))
        continue
    end
    %}
    
    
    if(isempty(currThetaSeqData.thetaSequenceInfoPerSequence))
        continue
    end
    
    %{
    if(~isfield(currThetaSeqData,'thetaSequenceInfoPerSequenceAtLeast4Fields'))
        continue
    end
    %}
    
    thetaSequenceInfoPerSequence=currThetaSeqData.thetaSequenceInfoPerSequence;
       % thetaSequenceInfoPerSequence=currThetaSeqData.thetaSequenceInfoPerSequenceAtLeast4Fields;

    
    seqStrs=fieldnames(thetaSequenceInfoPerSequence);
    
    for ti=1:length(seqStrs)
        
        if(exampleOnly && ~strcmp(seqStrs{ti},exampleFieldSeqStr))
            continue
        end
        currSeqData=thetaSequenceInfoPerSequence.(seqStrs{ti});
        numOccurences=currSeqData.numOccurences;
       
        avgFieldWidth=currSeqData.avgFieldWidth;
        numFieldsInSeq=length(currSeqData.fieldCenterSeq);
        travSpeedPerOcc=cell2num(currSeqData.travAvgSpeedPerOcc);
        
        distColors=jet(numFieldsInSeq);
                 
                 if(numFieldsInSeq==2)
                     distColors=[ 0     0     1; 1     0     0];
                 end
                 
                 if(numFieldsInSeq==3)
                     distColors=[ 0     0     1;0     1     0; 1     0     0];
                 end
        
        seqDirStr=currSeqData.seqDirStr;
        
        [sortedSpeedPerOcc,sortBySpeedOccIdxes]=sort(travSpeedPerOcc);
        
      
        allPhaseSeqPerSpeedSortedOcc=currSeqData.allPhaseSeqPerOcc(sortBySpeedOccIdxes);
        meanPhaseSeqPerSpeedSortedOcc=currSeqData.meanPhaseSeqPerOcc(sortBySpeedOccIdxes);
        
        allPhaseSeqPerOriginalOcc=currSeqData.allPhaseSeqPerOcc;
         meanPhaseSeqPerOriginalOcc=currSeqData.meanPhaseSeqPerOcc;
         
        timeInFieldSeqPerOriginalOcc=currSeqData.timeSinceFieldStartSeqPerOcc;
        timeToFieldEndSeqPerOriginalOcc=currSeqData.timeToFieldEndSeqPerOcc; %faster=less time remaining
        
        fieldDurationSeqPerOriginalOcc=currSeqData.fieldDurationSeqPerOcc;
        
        timeOfFieldInSeqPerOccPerFieldInSec=currSeqData.timeOfFieldInSeqPerOccPerFieldInSec;
        
        timeSinceTravStartEveryCyclePerOccInSec=currSeqData.timeSinceTravStartEveryCyclePerOcc;
        
        totalTraversalDurationPerOcc=currSeqData.totalTraversalDurationPerOcc;
        
        timeOfFieldInSeqPerOccPerFieldInSec=[timeOfFieldInSeqPerOccPerFieldInSec{:,:}];
        
        fieldDurationSeqPerFieldPerOcc=[fieldDurationSeqPerOriginalOcc{:,:}];
        timeSinceTravStartEveryCyclePerOccInSec=[timeSinceTravStartEveryCyclePerOccInSec{:,:}];
        
          gaussKernelSDduration=0.3; 
        durDistrEdges=linspace(0,10,100);
        
        totalTraversalDurationPerOcc=[totalTraversalDurationPerOcc{:,:}];
        
           %fH=figure;
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           %typical duration as peak duration
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       [currPdist,typicalTotalTraversalDuration] = getKernelDistrAndPeakVal(totalTraversalDurationPerOcc,gaussKernelSDduration,durDistrEdges);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           %typical duration as mean duration
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %typicalTotalTraversalDuration=nanmean(totalTraversalDurationPerOcc);
        
        timeSinceTravStartEveryCyclePerOccSeqFrac=timeSinceTravStartEveryCyclePerOccInSec/typicalTotalTraversalDuration;
        
       %timeOfFieldInSeqPerOccPerFieldFrac=timeOfFieldInSeqPerOccPerFieldInSec/typicalTotalTraversalDuration;
       
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           %normalize to slowest travesral
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       maxTimeOfFieldInSeq=max(timeOfFieldInSeqPerOccPerFieldInSec(:));
       
              timeOfFieldInSeqPerOccPerFieldFrac=timeOfFieldInSeqPerOccPerFieldInSec/maxTimeOfFieldInSeq;
                            %timeOfFieldInSeqPerOccPerFieldFrac=timeOfFieldInSeqPerOccPerFieldInSec/typicalTotalTraversalDuration;

             % typicalTotalTraversalDuration

        typicalFieldDurationPerField=NaN(numFieldsInSeq,1);
     
        %gaussKernelSDduration=1;
        
        normalizeByLongestDuration=1;
         normalizeByLongestDuration=0;
        
        for fdi=1:numFieldsInSeq
            currFieldAllDurations=fieldDurationSeqPerFieldPerOcc(fdi,:);
            
            [currPdist,peakVal] = getKernelDistrAndPeakVal(currFieldAllDurations,gaussKernelSDduration,durDistrEdges);
            if(normalizeByLongestDuration)
                typicalFieldDurationPerField(fdi)=max(currFieldAllDurations);
            else
                typicalFieldDurationPerField(fdi)=peakVal;
            end
        end
         %figure; plot(fieldDurationSeqPerOccPerFieldArray(1,:),fieldDurationSeqPerOccPerFieldArray(2,:),'k.','MarkerSize',20)
        
        fieldMaxRatePosPerField=NaN(numFieldsInSeq,1);

        for ai=1:numFieldsInSeq
            currFieldFiringRatePerPos=currSeqData.firingRatePerPosPerFieldInSeq{ai};

            [~,currFieldMaxRatePosBin]=max(smooth(currFieldFiringRatePerPos(:,2)));

            currFieldMaxRatePos=currFieldFiringRatePerPos(currFieldMaxRatePosBin,1);

            fieldMaxRatePosPerField(ai)=currFieldMaxRatePos;
        end
        
        if(strcmp(seqDirStr,'rightward'))
            [~,fieldIdxToCenterRank]=sort(fieldMaxRatePosPerField);
        else(strcmp(seqDirStr,'leftward'))
            [~,fieldIdxToCenterRank]=sort(fieldMaxRatePosPerField,'descend');
        end
                
        if(numFieldsInSeq<minNumFieldsInSeq)
            continue
        end
        
        
        
        speedStdDevAcrossLaps=nanstd(cell2num(currSeqData.travAvgSpeedPerOcc));
        
        if(numOccurences>=minNumOccurences && speedStdDevAcrossLaps>=minSpeedStdDevAcrossLaps...
                & numFieldsInSeq>=minNumFieldsInSeq)
            %currSeqData
            %fH=figure; 
        distFracColors=jet(numFieldsInSeq);
        allFieldDistFracs=[currSeqData.fieldDistFracSeqPerOcc{:}];
        
        %repFracTol=0.05;
        %repFracTol=0.2/3;
        %repFracTol=0.075;
        %repFracTol=0.025;
        %repFracTol=0.1;
        %repFracTol=0.075;
        repFracTol=1/12; %total 1/6th of field dist tolerance
        %repFracTol=0.025;
        %repFracTol=0.1;
         %repFracTol=0.15;
         
        isControlledDistPerOcc=ones(numOccurences,1);
        
        minEdgeDistFrac=0.1;
        maxEdgeDistFrac=0.9;
        for fi=1:numFieldsInSeq
            currFieldDistFracs=allFieldDistFracs(fi,:);
            
            gaussKernelSD=0.2/3; %3 of these is approximately full width
            %gaussKernelSD=0.15/3; %3 of these is approximately full width
            %gaussKernelSD=0.3/3; %3 of these is approximately full width
            
           
            if(noEdges)
                distrEdges=linspace(minEdgeDistFrac,maxEdgeDistFrac,100);
            else
                distrEdges=linspace(0,1,100);
            end
            
            [currPdist,currFieldRepFrac] = getKernelDistrAndPeakVal(currFieldDistFracs,gaussKernelSD,distrEdges);
            
            currFieldCloseToRepFracOccIdxes=currFieldDistFracs>=currFieldRepFrac-repFracTol & currFieldDistFracs<=currFieldRepFrac+repFracTol;
            
            isControlledDistPerOcc=isControlledDistPerOcc(:) & (currFieldCloseToRepFracOccIdxes(:)); %must be good for all fields
            
            if(noEdges)
                isControlledDistPerOcc(currFieldDistFracs<minEdgeDistFrac | currFieldDistFracs>maxEdgeDistFrac)=0;
            end
            
            currFieldDistFracs(~isControlledDistPerOcc)=NaN;
            
            %{
            hold on
            plot(1:length(currFieldDistFracs),currFieldDistFracs,'.','Color',distFracColors(fi,:),'MarkerSize',30)
             hold on
             hline(currFieldRepFrac,'k--',3)
            %}
            
            
        end
        
        if(~controlForDistanceAcrossLaps)
             isControlledDistPerOcc=ones(numOccurences,1);
        end
        meanDistFracPerField=nanmean([currSeqData.fieldDistFracSeqPerOcc{:}],2);
       
        %isControlledDistPerOcc
        %hline(meanDistFracPerField,'k-',3)
        
        %close all
        %continue
        
        
        
        if(usePeakPhaseAsSummary)
            peakPhaseSeqPerOcc=currSeqData.peakPhaseSeqPerOcc;
        else
            
              peakPhaseSeqPerOcc=currSeqData.meanPhaseSeqPerOcc;
        end

            %travAvgSpeedPerOcc=currSeqData.travAvgSpeedPerOcc;
            %travSpeedPerOcc=cell2num(travAvgSpeedPerOcc);
            
            if(speedInFieldsPerSec)
                travSpeedPerOcc=travSpeedPerOcc/avgFieldWidth; %from m/s to field widths per second
            end
            
            if(speedInFieldsPerSec)
                sortedSpeedPerOcc=sortedSpeedPerOcc/avgFieldWidth;
            end
            travSpeedTracePerOcc=currSeqData.travSpeedTracePerOcc;
            numFieldsInSeq=length(currSeqData.fieldCenterSeq);
            
            avgPhaseDiffPerOcc=NaN(numOccurences,1);
            avgPhasePerOcc=NaN(numOccurences,1);
            speedVarWithinTrav=NaN(numOccurences,1);
            
            for oi=1:numOccurences
                
                if(~isControlledDistPerOcc(oi))
                    continue
                end
                %for fi=1:numFieldsInSeq
                    currOccPhaseSeq=peakPhaseSeqPerOcc{oi};
                    
                    %avgPhaseDiffPerOcc(oi)=nanmean(angdiffDeg(currOccPhaseSeq)); %measure of phase spread
                   if(useCircMeasureSpread)
                         avgPhaseDiffPerOcc(oi)=1-circ_r(ang2rad(currOccPhaseSeq(:)));
                   else
                         avgPhaseDiffPerOcc(oi)=circMeanDegNoMod360(angdiffDeg(currOccPhaseSeq)); %measure of phase spread

                   end

                    if(nonCircStats)
                        avgPhasePerOcc(oi)=nanmean(currOccPhaseSeq);             %measure of central phase
                    else
                        avgPhasePerOcc(oi)=circMeanDeg(currOccPhaseSeq);             %measure of central phase
                    end

                     %avgPhasePerOcc(oi)=circMedianDeg(currOccPhaseSeq);             %measure of central phase
                     
                     speedVarWithinTrav(oi)=nanstd(travSpeedTracePerOcc{oi});
                     currTravSpeedTrace=travSpeedTracePerOcc{oi};
                     nonNanSpeedTrace=currTravSpeedTrace(~isnan(currTravSpeedTrace));
                     %speedVarWithinTrav(oi)=getSEMacrossRows(nonNanSpeedTrace(:));
                     
                     
                %end
            end
            
            
            
            computeCorrIdxes=travSpeedPerOcc>=minSpeed & travSpeedPerOcc<=maxSpeed;
            
            constSpeedTravIdxes=speedVarWithinTrav <= maxSpeedVarWithinTrav;
            
            %Dec22
            %{
            avgPhaseDiffPerOcc(~constSpeedTravIdxes)=NaN;
            avgPhasePerOcc(~constSpeedTravIdxes)=NaN;
            %}
            
         
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %absolute cutoff
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            highSpeedMinCutoff=0.5;
            lowSpeedMaxCutoff=0.3;
            
             highSpeedMinCutoff=0.45;
            lowSpeedMaxCutoff=0.35;
            
            highSpeedMinCutoff=0.55;
            lowSpeedMaxCutoff=0.25;
            
             % highSpeedMinCutoff=0.4;
            %lowSpeedMaxCutoff=0.4;
            
             %highSpeedMinCutoff=0.5;
            %lowSpeedMaxCutoff=0.3;
            
             highSpeedMinCutoff=0.45;
            lowSpeedMaxCutoff=0.35;
             %highSpeedMinCutoff=0.6;
            %lowSpeedMaxCutoff=0.2;
             highSpeedMinCutoff=0.5;
            lowSpeedMaxCutoff=0.3;
            
            highSpeedMinCutoff=0.55;
            lowSpeedMaxCutoff=0.25;
            
            if(speedInFieldsPerSec) 
                highSpeedMinCutoff=1.1;
                %lowSpeedMaxCutoff=0.5;
                %highSpeedMinCutoff=0.8;
                %lowSpeedMaxCutoff=0.7;
                %lowSpeedMaxCutoff=0.8;
                lowSpeedMaxCutoff=0.75;
                
                highSpeedMinCutoff=1;
                lowSpeedMaxCutoff=0.4;
                
                 highSpeedMinCutoff=1;
                lowSpeedMaxCutoff=0.5;
                
                highSpeedMinCutoff=1.1; %p ***
                lowSpeedMaxCutoff=0.55;
                
                 highSpeedMinCutoff=1.05; %p **
                lowSpeedMaxCutoff=0.55;
                
                %highSpeedMinCutoff=1.05; %p ***
                %lowSpeedMaxCutoff=0.6;
            end
            
            %highSpeedMinCutoff=0.6;
            %lowSpeedMaxCutoff=0.2;
            %highSpeedIdxes=computeCorrIdxes & travSpeedPerOcc>=highSpeedMinCutoff;
            %lowSpeedIdxes=computeCorrIdxes & travSpeedPerOcc<=lowSpeedMaxCutoff;
            
             %highSpeedIdxes=travSpeedPerOcc>=highSpeedMinCutoff;
            %lowSpeedIdxes=travSpeedPerOcc<=lowSpeedMaxCutoff;
            
            %top and bottom 33 percentile as high speed and low speed
            %highSpeedMinCutoff=prctile(travSpeedPerOcc,2/3*100);
            %lowSpeedMaxCutoff=prctile(travSpeedPerOcc,1/3*100);
            
            %highSpeedMinCutoff=prctile(travSpeedPerOcc,80);
            %lowSpeedMaxCutoff=prctile(travSpeedPerOcc,20);
            
            %highSpeedMinCutoff=prctile(travSpeedPerOcc,90);
            %lowSpeedMaxCutoff=prctile(travSpeedPerOcc,10);
            
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %relative to ind cell speed cutoff
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            highSpeedMinCutoff=prctile(travSpeedPerOcc,indFieldHighSpeedPrctileCutoff);
            lowSpeedMaxCutoff=prctile(travSpeedPerOcc,indFieldLowSpeedPrctileCutoff);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %absolute speed limits to ensure sufficient variation within
            %individual field
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             % highSpeedMinCutoff=0.7;
             %lowSpeedMaxCutoff=0.5;
             % highSpeedMinCutoff=0;
             %lowSpeedMaxCutoff=Inf;
            
            highSpeedIdxes=constSpeedTravIdxes(:) & travSpeedPerOcc(:)>=highSpeedMinCutoff;
            lowSpeedIdxes=constSpeedTravIdxes(:) & travSpeedPerOcc(:)<=lowSpeedMaxCutoff;
            midSpeedIdxes=constSpeedTravIdxes(:) & travSpeedPerOcc(:)>lowSpeedMaxCutoff & travSpeedPerOcc(:)<highSpeedMinCutoff;
            
            highSpeedDistFracs=currFieldDistFracs(highSpeedIdxes);
            lowSpeedDistFracs=currFieldDistFracs(lowSpeedIdxes);
            
            %figure place holder
            circFigH=[];
            
            %enoughPtsInEachCategory=sum(highSpeedIdxes)>=minNumPtsInSpeedCategory && sum(lowSpeedIdxes)>=minNumPtsInSpeedCategory;
                  
            enoughPtsInEachCategory=sum(highSpeedIdxes)+sum(lowSpeedIdxes)>=minNumPtsInSpeedCategory*2;
            
            if(~enoughPtsInEachCategory)
                %continue
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %plot summary of this field sequence vs speed
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %if(showPlots && enoughPtsInEachCategory)
            showSpeedSorted2DthetaSeq=1;
            if(showPlots)
                
                figure
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 %rate fields
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                subplot(3,3,1)
                %distColors=jet(numFieldsInSeq);

                
                for ai=1:numFieldsInSeq
                    currFieldFiringRatePerPos=currSeqData.firingRatePerPosPerFieldInSeq{ai};
                    plot(currFieldFiringRatePerPos(:,1),currFieldFiringRatePerPos(:,2),'Color',distColors(fieldIdxToCenterRank(ai),:),'LineWidth',3)
                    hold on
                end
                
                box off
                
                %subplot(3,3,1)
                
                for oi=1:numOccurences
                            hold on
                            
                            if(~isControlledDistPerOcc(oi)) %skip if not in controlled distance range
                               continue
                            end
                           
                            if(strcmp(seqDirStr,'rightward'))
                                vline(currSeqData.positionPerOcc{oi},'k--',1)
                            else
                                maxTrackLength=2.5;
                                vline(maxTrackLength-currSeqData.positionPerOcc{oi},'k--',1)
                            end
                end
                            
                
                xlabel('Position (m)')
              ylabel('Firing rate (Hz)')
              axis tight
                
                %colormap(jet(numFieldsInSeq))
                colormap(distColors)
               cb=colorbar;
               cb.Ticks=[];
                ylabel(cb,'Behavioral field order')
                caxis([1 numFieldsInSeq])
                
                title(sprintf('Place fields in behavioral sequence'))%, %s %s',currSessName,removeUnderscores(seqStrs{ti})))
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 %distance vs speed
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                subplot(3,3,2)
 
                 

                
                for ai=1:numFieldsInSeq
                    currFieldDistFracs=allFieldDistFracs(ai,:);
                    if(controlForDistanceAcrossLaps)
                        plot(travSpeedPerOcc(isControlledDistPerOcc),currFieldDistFracs(isControlledDistPerOcc),'.', 'Color',distColors(fieldIdxToCenterRank(ai),:),'MarkerSize',30)
                    else
                        plot(travSpeedPerOcc,currFieldDistFracs,'.', 'Color',distColors(fieldIdxToCenterRank(ai),:),'MarkerSize',30)
                    end
                    hold on
                end
                
                box off
                
                xlim([minSpeed maxSpeed])
                ylim([0 1])
                xlabel('Running speed (field/sec)')
              ylabel('Distance in field (frac)')
              %axis tight
                
                %colormap(jet(numFieldsInSeq))
                colormap(distColors)
               cb=colorbar;
               cb.Ticks=[];
                ylabel(cb,'Behavioral field order')
                caxis([1 numFieldsInSeq])
                
                title(sprintf('Distance in field vs speed'))%, %s %s',currSessName,removeUnderscores(seqStrs{ti})))
               
                
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 %speed variability
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                subplot(3,3,3)
 
                %speedEdges=linspace(minSpeed,maxSpeed,21);
                speedEdges=linspace(minSpeed,maxSpeed,31);
                histogram(travSpeedPerOcc(isControlledDistPerOcc),speedEdges)
                hold on

                 pLh=vline(lowSpeedMaxCutoff,'k--',3);
                    pHh=vline(highSpeedMinCutoff,'k-',3);
                    
                    %legend([pLh pHh],{sprintf('%.2f fields/sec',lowSpeedMaxCutoff),sprintf('%.2f fields/sec',highSpeedMinCutoff)},'Location','westoutside')
                   %legend([pLh pHh],{sprintf('%.2f fields/sec',lowSpeedMaxCutoff),sprintf('%.2f fields/sec',highSpeedMinCutoff)},'Location','northeast')
                   legend([pLh pHh],{sprintf('%.2f (%dth prctile)',lowSpeedMaxCutoff, round(indFieldLowSpeedPrctileCutoff)),sprintf('%.2f (%dth prctile)',highSpeedMinCutoff,round(indFieldHighSpeedPrctileCutoff))},'Location','best')


                    legend boxoff

                 %pS=vline(currMeanTravSpeed,'k--',2);
             
             %legend([pS],{'curr. speed'},'Location','northeast')%'eastoutside')
             %legend boxoff
            
             xlim([minSpeed,maxSpeed])
             
             if(speedInFieldsPerSec)
                 xlabel('traversal speed (fields/sec)')
             else
                xlabel('traversal speed (m/s)')
             end
            
            ylabel('traversal count')
             axis tight
             %axis square
             box off

             title('Speed distribution across laps')
             hold off
           end
             
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %speed sorted 2D theta sequences
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if(showSpeedSorted2DthetaSeq)% && enoughPtsInEachCategory)
                    
                    phaseOnXaxis=1;
                    phaseOnXaxis=0;
                    
                    %subplot(1,2,2)
                    %subplot(2,3,[4 5 6])
                    %subplot(2,3,[3])
                    highSpeedOccCount=0;
                    lowSpeedOccCount=0;
                    midSpeedOccCount=0;

                    %use these if set tight subplots
                     %axisHighSpeed=subplot(3,2,3);
                     %axisLowSpeed=subplot(3,2,4);
                     
                     lowSpeedFieldFracPerFieldPerOcc=NaN(numFieldsInSeq,numOccurences);
                     lowSpeedMeanPhasePerFieldPerOcc=NaN(numFieldsInSeq,numOccurences);
                     
                     midSpeedFieldFracPerFieldPerOcc=NaN(numFieldsInSeq,numOccurences);
                     midSpeedMeanPhasePerFieldPerOcc=NaN(numFieldsInSeq,numOccurences);
                       
                     highSpeedFieldFracPerFieldPerOcc=NaN(numFieldsInSeq,numOccurences);
                     highSpeedMeanPhasePerFieldPerOcc=NaN(numFieldsInSeq,numOccurences);
                     
                     timeOfFieldToSeqEndPerOccPerFieldFrac=1-timeOfFieldInSeqPerOccPerFieldFrac;
                     
                  %thetaSequenceInfoPerSequence.(seqStrs{ti}).timeOfFieldToSeqEndPerOccPerFieldFrac=timeOfFieldToSeqEndPerOccPerFieldFrac;
                  thetaSequenceInfoPerSequence.(seqStrs{ti}).timeOfFieldInSeqPerOccPerFieldFrac=timeOfFieldInSeqPerOccPerFieldFrac;


                  
                  thetaCycleWaveformPerOcc=currSeqData.thetaCycleWaveform;
                     
                    for oi=1:numOccurences
                        
                        if(~isControlledDistPerOcc(oi))
                            continue
                        end
                            

                        
                        if(highSpeedIdxes(oi))
                            if(showPlots)
                                subplot(3,3,6)
                            end
                           %axes(axisHighSpeed)
                     
                            highSpeedOccCount=highSpeedOccCount+1;
                        elseif(lowSpeedIdxes(oi))
                            if(showPlots)
                                subplot(3,3,4)  
                            end
                            %axes(axisLowSpeed)
                            lowSpeedOccCount=lowSpeedOccCount+1;
                        elseif(midSpeedIdxes(oi))
                            if(showPlots)
                                subplot(3,3,5)  
                            end
                            %axes(axisLowSpeed)
                            midSpeedOccCount=midSpeedOccCount+1;
                        else
                            continue %not valid occurence (e.g. non-constant speed, uncontrolled distance)
                        end
                        
                        currOccAllPhaseSeq=allPhaseSeqPerOriginalOcc{oi};
                        currOccTimeInFieldSeq=timeInFieldSeqPerOriginalOcc{oi};
                        currOccTimeToFieldEndSeq=timeToFieldEndSeqPerOriginalOcc{oi};
                        currOccFieldDurationSeq=fieldDurationSeqPerOriginalOcc{oi};
                        
                        currOccFieldDistFracSeq=allFieldDistFracs(:,oi);
                        
                        %currOccTimeOfFieldToSeqEndFracPerField=1-timeOfFieldInSeqPerOccPerFieldFrac(:,oi);
                        currOccTimeOfFieldToSeqEndFracPerField=timeOfFieldToSeqEndPerOccPerFieldFrac(:,oi);
                        currOccFieldDurationPerField=fieldDurationSeqPerOriginalOcc{oi};
                         for ofi=1:length(currOccAllPhaseSeq)
                                currOccFieldPhases=currOccAllPhaseSeq{ofi};
                                %currOccFieldFracs=currOccTimeInFieldSeq(ofi)/currOccFieldDurationSeq(ofi);
                                %currOccFieldFracs=currOccTimeOfFieldToSeqEndFracPerField(ofi);
                                %currOccFieldFracs=currOccTimeToFieldEndSeq(ofi)/typicalFieldDurationPerField(ofi);
                                %currOccFieldFracs=(currOccTimeInFieldSeq(ofi)/typicalFieldDurationPerField(ofi));
                                 %currOccFieldFracs=1-(currOccTimeToFieldEndSeq(ofi)/nanmean(typicalFieldDurationPerField));
                                %currOccFieldFracs=(currOccTimeToFieldEndSeq(ofi)/nanmean(typicalFieldDurationPerField));
                                 %currOccFieldFracs=(currOccTimeToFieldEndSeq(ofi)/max(typicalFieldDurationPerField));
                                 %currOccFieldFracs=1-(currOccTimeInFieldSeq(ofi)/max(typicalFieldDurationPerField));
                                 %currOccFieldFracs=1-(currOccTimeInFieldSeq(ofi)/currOccFieldDurationPerField(ofi));
                                 %currOccFieldFracs=-(currOccTimeInFieldSeq(ofi)/max(typicalFieldDurationPerField));

                                 %time to end position, if maintaining avg speed from current position
                                 %remaining distance over speed, relative
                                 %to typical duration
                                 %currOccFieldFracs=((1-currOccFieldDistFracSeq(ofi))/travSpeedPerOcc(oi))/max(typicalFieldDurationPerField) ;
                                 currOccFieldFracs=(currOccTimeToFieldEndSeq(ofi))/max(typicalFieldDurationPerField) ;

                                                                  
                                 currOccThetaWaveformPerTime=thetaCycleWaveformPerOcc{oi};
                                 currOccThetaWaveform=currOccThetaWaveformPerTime(:,2);
                                 
                                 thetaWaveformPhaseAxis=linspace(0,360,length(currOccThetaWaveform));
                                

                                %plotRasterStyle(currOccFieldPhases,lowSpeedOccCount,NaN,NaN,distColors(fieldIdxToCenterRank(ofi),:))
                                
                                if(showPlots)
                                    if(phaseOnXaxis)
                                        plot(currOccFieldPhases,currOccFieldFracs,'.','Color',distColors(fieldIdxToCenterRank(ofi),:),'MarkerSize',15)
                                    else
                                        plot(currOccFieldFracs,currOccFieldPhases,'.','Color',distColors(fieldIdxToCenterRank(ofi),:),'MarkerSize',15)
                                    end

                                    hold on


                                    %rescaledThetaWaveform=scaledata(currOccThetaWaveform,-1,0);
                                    rescaledThetaWaveform=currOccThetaWaveform;


                                    if(showThetaWaveforms)
                                        yyaxis right
                                        plot(thetaWaveformPhaseAxis,rescaledThetaWaveform,'k-')
                                        ylim([-2.5 2.5])
                                        ylabel('LFP (Z)')
                                        yyaxis left
                                    end

                                end
                                 %xDispUpLim=1.2;
                                %xDispLowLim=-0.2;
                                 xDispUpLim=1;
                                xDispLowLim=0;
                                %xDispUpLim=0;
                                %xDispLowLim=-1;
                                
                                %validIdxes=currOccFieldFracs>=0 & currOccFieldFracs<=1;
                                %validIdxes=currOccFieldFracs>-Inf & currOccFieldFracs<Inf;
                                validIdxes=currOccFieldFracs>=xDispLowLim & currOccFieldFracs<=xDispUpLim;
                                
                                if(highSpeedIdxes(oi))
                                    highSpeedFieldFracPerFieldPerOcc(ofi,oi)=nanmean(currOccFieldFracs(validIdxes));
                                    highSpeedMeanPhasePerFieldPerOcc(ofi,oi)=circMeanDeg(currOccFieldPhases(validIdxes));
                                elseif(lowSpeedIdxes(oi))
                                    lowSpeedFieldFracPerFieldPerOcc(ofi,oi)=nanmean(currOccFieldFracs(validIdxes));
                                    lowSpeedMeanPhasePerFieldPerOcc(ofi,oi)=circMeanDeg(currOccFieldPhases(validIdxes));
                                elseif(midSpeedIdxes(oi))
                                    midSpeedFieldFracPerFieldPerOcc(ofi,oi)=nanmean(currOccFieldFracs(validIdxes));
                                    midSpeedMeanPhasePerFieldPerOcc(ofi,oi)=circMeanDeg(currOccFieldPhases(validIdxes));
                                end
                                
                                
                       
                                hold on
                                
                         end

                         if(showPlots)
                             
                            if(highSpeedIdxes(oi))
                                 title('high speed traversal 2D theta sequences')
                            elseif(lowSpeedIdxes(oi))
                                 title('low speed traversal 2D theta sequences')
                            elseif(midSpeedIdxes(oi))
                                title('mid speed traversal 2D theta sequences')
                            end

                                  box off

                              timeInFieldDisp=linspace(0,1,100);
                              [logTimeToEndDisp] = getLogTime(timeInFieldDisp,1.4);
                              %[logTimeToEndDisp] = getLogTime(1-timeInFieldDisp,1.4);
                              %plot(timeInFieldDisp,logTimeToEndDisp*360,'k--','LineWidth',3)
                              %plot(1-timeInFieldDisp,(logTimeToEndDisp)*360,'k--','LineWidth',3)
                              %plot(timeInFieldDisp,(logTimeToEndDisp)*360,'r--','LineWidth',3)

                                if(phaseOnXaxis)
                                    plot((logTimeToEndDisp)*360,timeInFieldDisp,'k--','LineWidth',3)
                                else
                                    plot(timeInFieldDisp,(logTimeToEndDisp)*360,'k--','LineWidth',3)
                                end


                              %xDispUpLim=1.2;
                              %xDispLowLim=-0.2;

                              %xDispUpLim=1.5;
                              %xDispLowLim=-0.5;
                              %xDispUpLim=1;
                              %xDispLowLim=0;

                              %xDispUpLim=1;
                              %xDispLowLim=0;

                              if(phaseOnXaxis)
                                   ylim([xDispLowLim xDispUpLim])
                                  xlim([0 360])
                                   %daspect([360 (xDispUpLim-xDispLowLim) 1])

                                %ylabel('Time to field start (norm.)')
                                ylabel('Time to field end (norm.)') %time to end position, if maintaining avg speed from current position
                                xlabel('Theta phase (deg)')
                              else
                                xlim([xDispLowLim xDispUpLim])
                                  ylim([0 360])
                                  %daspect([1.2 360 1])
                                  if(~showThetaWaveforms)
                                   daspect([xDispUpLim-xDispLowLim 360 1])
                                  end

                                xlabel('Time to field end (frac)')
                                ylabel('Theta phase (deg)')

                              end

                               %colormap(jet(numFieldsInSeq))
                               colormap(distColors)
                               cb=colorbar;
                               cb.Ticks=[];
                                ylabel(cb,'Behavioral field order')
                                caxis([1 numFieldsInSeq])
                            
                         end
                    end
                    

                    for sli=1:3
                            if(showPlots)
                                 subplot(3,3,3+sli)         
                           
                            end
                           if(sli==1)
                               currSpeedLevelFieldFracPerFieldPerOcc=lowSpeedFieldFracPerFieldPerOcc;
                               currSpeedLevelMeanPhasePerFieldPerOcc=lowSpeedMeanPhasePerFieldPerOcc;
                           elseif(sli==2)
                               currSpeedLevelFieldFracPerFieldPerOcc=midSpeedFieldFracPerFieldPerOcc;
                               currSpeedLevelMeanPhasePerFieldPerOcc=midSpeedMeanPhasePerFieldPerOcc;
                           elseif(sli==3)
                               currSpeedLevelFieldFracPerFieldPerOcc=highSpeedFieldFracPerFieldPerOcc;
                               currSpeedLevelMeanPhasePerFieldPerOcc=highSpeedMeanPhasePerFieldPerOcc;
                           end
                           
                        for mfi=1:numFieldsInSeq
                            currFieldAllOccFieldFrac=currSpeedLevelFieldFracPerFieldPerOcc(mfi,:);
                            currFieldAllOccPhases=currSpeedLevelMeanPhasePerFieldPerOcc(mfi,:);

                            currFieldMeanFieldFrac=nanmean(currFieldAllOccFieldFrac);
                            currFieldMeanPhase=circMeanDeg(currFieldAllOccPhases(~isnan(currFieldAllOccPhases)));

                            if(showPlots)
                                if(phaseOnXaxis)
                                    scatter(currFieldMeanPhase,currFieldMeanFieldFrac,100,'MarkerFaceColor',distColors(fieldIdxToCenterRank(mfi),:),'MarkerEdgeColor','k','LineWidth',1)
                                    hold on
                                    plot([0 currFieldMeanPhase],[currFieldMeanFieldFrac currFieldMeanFieldFrac],'-','Color',distColors(fieldIdxToCenterRank(mfi),:),'LineWidth',3);
                                    plot([currFieldMeanPhase currFieldMeanPhase],[currFieldMeanFieldFrac 0],'-','Color',distColors(fieldIdxToCenterRank(mfi),:),'LineWidth',3);

                                else
                                    scatter(currFieldMeanFieldFrac,currFieldMeanPhase,100,'MarkerFaceColor',distColors(fieldIdxToCenterRank(mfi),:),'MarkerEdgeColor','k','LineWidth',1)
                                    hold on
                                    plot([currFieldMeanFieldFrac currFieldMeanFieldFrac],[0 currFieldMeanPhase],'-','Color',distColors(fieldIdxToCenterRank(mfi),:),'LineWidth',3);
                                    %plot([currFieldMeanFieldFrac xDispUpLim],[currFieldMeanPhase currFieldMeanPhase],'-','Color',distColors(fieldIdxToCenterRank(mfi),:),'LineWidth',3);
                                    plot([0 currFieldMeanFieldFrac],[currFieldMeanPhase currFieldMeanPhase],'-','Color',distColors(fieldIdxToCenterRank(mfi),:),'LineWidth',3);
                                end
                            end

                        end
                    
                    end
                    
 
                       
                       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                       %individual field sequence speed vs phase signature plots
                       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                       
                       useModifiedCircStats=1;
                       useModifiedCircStats=0;

                    allUsedRunSpeedFieldsPerSec=NaN(numOccurences,1);
                    allUsedDistFracs=NaN(numOccurences,1);
                    allUsedDistSepFracs=NaN(numOccurences,1);
                    allUsedTimeSepFracs=NaN(numOccurences,1);
                    allUsedPhaseSeparations=NaN(numOccurences,1);
                    allUsedMeanPhases=NaN(numOccurences,1);
                    allUsedPhaseTimeCompressions=NaN(numOccurences,1);
                    
                    %usePhaseOffsetAsEarlyFieldPhase=1;
                     usePhaseOffsetAsEarlyFieldPhase=0;
                    
                     
                     
                       for oi=1:numOccurences
                           
                           if(~isControlledDistPerOcc(oi)) %skip if not in controlled distance range
                               continue
                           end
                           
                         
                           if(~constSpeedTravIdxes(oi))
                               continue
                           end
                            currOccMeanPhaseSeq=meanPhaseSeqPerOriginalOcc{oi};
                            currOccFieldDistFracSeq=allFieldDistFracs(:,oi);
  
                            currOccTimeInFieldSeqSec=timeInFieldSeqPerOriginalOcc{oi};
                            
                            %currOccTimeToFieldEndSeq=timeToFieldEndSeqPerOriginalOcc{oi}
                            
                            %currOccTimeToFieldEndSeqSec=timeToFieldEndSeqPerOriginalOcc{oi};
                            %currOccTimeToFieldEndSeq=1-(currOccTimeInFieldSeqSec/nanmean(typicalFieldDurationPerField));
                            %currOccTimeToFieldEndSeq=1-(currOccTimeInFieldSeqSec/max(typicalFieldDurationPerField));
                            
                            %currOccTimeToFieldEndSeq=1-(currOccTimeInFieldSeqSec./(fieldDurationSeqPerOriginalOcc{oi}));
                            %currOccTimeToFieldEndSeq=-currOccTimeInFieldSeqSec/max(typicalFieldDurationPerField);
                            
                            %currOccTimeToFieldEndSeq=((1-currOccFieldDistFracSeq)/travSpeedPerOcc(oi))/max(typicalFieldDurationPerField) ;

                            currOccTimeToFieldEndSeqFrac=currOccTimeInFieldSeqSec/max(typicalFieldDurationPerField) ;
                            
                            currOccRunSpeedFieldsPerSec=travSpeedPerOcc(oi);
                            
                            %SWITCH FROM CIRC DIFF TO LINEAR DIFF JAN 15 2022
                            currOccPhaseSeparationSeq=angdiffDeg(currOccMeanPhaseSeq);
                            %currOccPhaseSeparationSeq=diff(currOccMeanPhaseSeq);
                            
                            currOccPhaseSeparation=circMeanDegNoMod360(currOccPhaseSeparationSeq);
                            %currOccPhaseSeparation=nanmean(currOccPhaseSeparationSeq);

                            %currOccPhaseSeparation=nanmean(angdiffDeg(currOccMeanPhaseSeq));

                               
                            
                            if(useModifiedCircStats && currOccPhaseSeparation<0) %if separation is beyond 180 degrees, write in positive form and use non circ mean (treated as linear variable within theta sequence)
                                currOccPhaseSeparation=currOccPhaseSeparation+360;
                                currOccMeanPhase=nanmean(currOccMeanPhaseSeq);
                            else
                                 currOccMeanPhase=circMeanDeg(currOccMeanPhaseSeq);
                            end
                            %{
                             %SWITCH FROM CIRC DIFF TO LINEAR DIFF JAN 15
                             2022 - not appropriate
                            currOccMeanPhase=nanmean(currOccMeanPhaseSeq);
                            %}
                            
                            
                              if(usePhaseOffsetAsEarlyFieldPhase)
                                    currOccMeanPhase=currOccMeanPhaseSeq(1); %earliest field mean phase
                                    
                                    %{
                                    if(currOccTimeInFieldSeqSec(1)<currOccTimeInFieldSeqSec(2)) %fields out of order
                                        disp('')
                                        
                                    end
                                    %}
                              end
                            
                            %currOccTimeSepSeq=(diff(currOccTimeToFieldEndSeqFrac));
                            
                            currOccTimeSepSeq=-(diff(currOccTimeToFieldEndSeqFrac));
                            
                            currOccMeanDistFrac=nanmean(currOccFieldDistFracSeq);
                            %currOccDistFracRange=range(currOccFieldDistFracSeq);
                            currOccAvgDistSepFrac=abs(nanmean(diff(currOccFieldDistFracSeq)));
                            
                            %currOccPhaseSepPerTimeSep=currOccPhaseSeparation/currOccTimeSep;
                              
                            
                            %currOccPhaseSepPerTimeSep=abs(currOccPhaseSeparation)/abs(currOccTimeSep); %magnitude of compression only
                             %currOccPhaseSepPerTimeSep=abs(nanmean(currOccPhaseSeparationSeq(:)./currOccTimeSepSeq(:)));
                          currOccPhaseSepPerTimeSep=abs(nanmean(currOccPhaseSeparationSeq(:)./currOccTimeSepSeq(:)));

                             %currOccPhaseSepPerTimeSep(currOccPhaseSepPerTimeSep<0)=NaN;
                             
                            %currOccPhaseSepPerTimeSep=abs(currOccPhaseSeparation)/abs(currOccTimeSep); %magnitude of compression only

                              
                            allUsedRunSpeedFieldsPerSec(oi)=currOccRunSpeedFieldsPerSec;
                            allUsedDistFracs(oi)=currOccMeanDistFrac;
                            allUsedDistSepFracs(oi)=currOccAvgDistSepFrac;
                            allUsedTimeSepFracs(oi)=nanmean(currOccTimeSepSeq);
                            
                            allUsedPhaseSeparations(oi)=currOccPhaseSeparation;
                             allUsedMeanPhases(oi)=currOccMeanPhase;
                             allUsedPhaseTimeCompressions(oi)=currOccPhaseSepPerTimeSep;
                            
                            if(currOccPhaseSepPerTimeSep==0)
                                disp('')
                            end
                           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                           %phase separation vs runnning speed per occurence
                           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                           if(showPlots)
                               %ptSize=15;
                               %ptSize=20;
                               ptSize=30;
                               vLineWidth=1;

                                subplot(3,3,7)
                                plot(currOccRunSpeedFieldsPerSec,currOccPhaseSeparation,'k.','MarkerSize',ptSize)
                                hold on

                                xlim([minSpeed maxSpeed])
                                if(useModifiedCircStats)
                                    ylim([0 360])
                                else
                                    ylim([-180 180])
                                end
                                hold on

                        pLh=vline(lowSpeedMaxCutoff,'k--',vLineWidth);
                        pHh=vline(highSpeedMinCutoff,'k-',vLineWidth);

                                xlabel('Running speed (fields/sec)')
                                ylabel('Phase separation (deg)')

                                title('Phase separation vs running speed')


                                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                               %mean phase vs runnning speed per occurence 
                               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                subplot(3,3,9)
                                plot(currOccRunSpeedFieldsPerSec,currOccMeanPhase,'k.','MarkerSize',ptSize)

                                hold on
                                xlim([minSpeed maxSpeed])
                                ylim([0 360])
                                xlabel('Running speed (fields/sec)')
                                ylabel('Mean phase (deg)')

                                title('Mean phase vs running speed')
                                hold on

                             pLh=vline(lowSpeedMaxCutoff,'k--',vLineWidth);
                                pHh=vline(highSpeedMinCutoff,'k-',vLineWidth);

                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                               %\deltaPhase/\deltaTime vs runnning speed per occurence
                               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                subplot(3,3,8)
                                plot(currOccRunSpeedFieldsPerSec,currOccPhaseSepPerTimeSep,'k.','MarkerSize',ptSize)


                                hold on

                                xlim([minSpeed maxSpeed])
                                %maxCompressionVal=1800;%up to five fields in a cycle
                                 maxCompressionVal=1080;%up to 3 field durations in a cycle

                                ylim([0 maxCompressionVal]) 
                                xlabel('Running speed (fields/sec)')
                                ylabel('phase-time compression (deg/field)')

                                hold on

                         pLh=vline(lowSpeedMaxCutoff,'k--',vLineWidth);
                        pHh=vline(highSpeedMinCutoff,'k-',vLineWidth);

                                title('Phase-time compression vs running speed')
                           end
                        end
                        
                        %hold off
                        
                        
                        numRunningAvgBins=5;
                        runningAvgAlpha=0.4;
                        speedEdgeBuffer=0.05;
                        runningAvgSpeedEdges=linspace(min(allUsedRunSpeedFieldsPerSec)-speedEdgeBuffer,max(allUsedRunSpeedFieldsPerSec)+speedEdgeBuffer,numRunningAvgBins+1);
                     [runningAvgSpeedBinCentersSep,phaseDiffRunningAvgValues] = getBinnedCircAverages(allUsedRunSpeedFieldsPerSec,allUsedPhaseSeparations,runningAvgSpeedEdges);
                     [runningAvgSpeedBinCentersMean,meanPhaseRunningAvgValues] = getBinnedCircAverages(allUsedRunSpeedFieldsPerSec,allUsedMeanPhases,runningAvgSpeedEdges);
                     [runningAvgSpeedBinCentersComp,phaseTimeCompRunningAvgValues] = getBinnedAverages(allUsedRunSpeedFieldsPerSec,allUsedPhaseTimeCompressions,runningAvgSpeedEdges);
                     
                     currFieldSeqRunSpeeds=allUsedRunSpeedFieldsPerSec(~isnan(allUsedRunSpeedFieldsPerSec));
                     currFieldSeqPhaseSep=allUsedPhaseSeparations(~isnan(allUsedPhaseSeparations));
                      currFieldSeqMeanDistFracs=allUsedDistFracs(~isnan(allUsedDistFracs));
                      currFieldSeqDistSepFracs=allUsedDistSepFracs(~isnan(allUsedDistSepFracs));
                      currFieldSeqTimeSepFracs=allUsedTimeSepFracs(~isnan(allUsedTimeSepFracs));
                      
                     currFieldSeqMeanPhase=allUsedMeanPhases(~isnan(allUsedMeanPhases));
                     currFieldSeqPhaseComp=allUsedPhaseTimeCompressions(~isnan(allUsedPhaseTimeCompressions));
                     %currFieldSeqPhaseComp=allUsedPhaseTimeCompressions(~isnan(allUsedMeanPhases));
                     allSessionRunSpeed=[allSessionRunSpeed;currFieldSeqRunSpeeds(:)];
                     allSessionMeanDistFracs=[allSessionMeanDistFracs; currFieldSeqMeanDistFracs(:)];
                     allSessionAvgDistSepFracs=[allSessionAvgDistSepFracs; currFieldSeqDistSepFracs(:)];
                     allSessionAvgTimeSepFracs=[allSessionAvgTimeSepFracs; currFieldSeqTimeSepFracs(:)];
                     
                     allSessionPhaseSep=[allSessionPhaseSep; currFieldSeqPhaseSep(:)];
                     allSessionMeanPhase=[allSessionMeanPhase; currFieldSeqMeanPhase(:)];
                     allSessionPhaseComp=[allSessionPhaseComp; currFieldSeqPhaseComp(:)];
                     
                     fieldSeqCount=fieldSeqCount+1;
                     
                     runSpeedsPerFieldSequence{fieldSeqCount}=currFieldSeqRunSpeeds;
                     phaseSepPerFieldSequence{fieldSeqCount}=currFieldSeqPhaseSep;
                     meanPhasePerFieldSequence{fieldSeqCount}=currFieldSeqMeanPhase;
                     phaseCompPerFieldSequence{fieldSeqCount}=currFieldSeqPhaseComp;
                     
                     meanDistPerFieldSequence{fieldSeqCount}=currFieldSeqMeanDistFracs;
                     distSepPerFieldSequence{fieldSeqCount}=currFieldSeqDistSepFracs;
                     timeSepPerFieldSequence{fieldSeqCount}=currFieldSeqTimeSepFracs;
                     
                     %meanDistFracPerFieldSequence{si}=currSessionMeanDistFracs;
                     %DistRangeFracPerFieldSequence{si}=currSessionDistRangeFracs;
                     
                     currMinSpeed=min(runningAvgSpeedBinCentersComp);
                                 currMaxSpeed=max(runningAvgSpeedBinCentersComp);
                                 currSpeedRange=range(runningAvgSpeedBinCentersComp);
                                 
                                 currMinDispSpeed=currMinSpeed-currSpeedRange*0.2;
                                 currMaxDispSpeed=currMaxSpeed+currSpeedRange*0.2;
                                 
                if(showPlots)
                         subplot(3,3,7)
                         box off
                         showRunningCircAvg=0;
                         %showRunningCircAvg=1;
                         showLinCircFit=1;

                             if(~isempty(phaseDiffRunningAvgValues))
                                 hold on
                                     if(showRunningCircAvg)
                                        pR=plot(runningAvgSpeedBinCentersSep,phaseDiffRunningAvgValues,'r-','LineWidth',2);
                                        pR.Color=[1 0 0 runningAvgAlpha];
                                     elseif(showLinCircFit)
                                         %[rho,p,slopeDegPerXunit, offsetDeg]=getCircCorrCoeff(runningAvgSpeedBinCentersSep,phaseDiffRunningAvgValues,fH);
                                        [rho,p,slopeDegPerXunit, offsetDeg]=getCircCorrCoeff( allUsedRunSpeedFieldsPerSec,allUsedPhaseSeparations,fH);

                                        
                                         allPhaseSepSpeedRho=[allPhaseSepSpeedRho;rho];
                                         allPhaseSepSpeedPval=[allPhaseSepSpeedPval;p];
                                     end
                                     %legend(pR,'running avg')

                             end

                          if(currMaxDispSpeed>currMinDispSpeed)
                                xlim([currMinDispSpeed currMaxDispSpeed])
                          end
                             
                             subplot(3,3,9)
                               box off
                             if(~isempty(meanPhaseRunningAvgValues))
                                 hold on
                                 if(showRunningCircAvg)
                                         pR=plot(runningAvgSpeedBinCentersMean,meanPhaseRunningAvgValues,'r-','LineWidth',2);

                                         %legend(pR,'running avg')
                                         pR.Color=[1 0 0 runningAvgAlpha];
                                 elseif(showLinCircFit)
                                         %[rho,p,slopeDegPerXunit, offsetDeg]=getCircCorrCoeff(runningAvgSpeedBinCentersSep,meanPhaseRunningAvgValues,fH);
                                         [rho,p,slopeDegPerXunit, offsetDeg]=getCircCorrCoeff(allUsedRunSpeedFieldsPerSec,allUsedMeanPhases,fH);

                                         
                                         allMeanPhaseSpeedRho=[allMeanPhaseSpeedRho;rho];
                                         allMeanPhaseSpeedPval=[allMeanPhaseSpeedPval;p];
                                 end


                             end
                          
                        if(currMaxDispSpeed>currMinDispSpeed)
                                xlim([currMinDispSpeed currMaxDispSpeed])
                          end
                         subplot(3,3,8)
                         if(~isempty(meanPhaseRunningAvgValues))
                              box off
                             hold on
                             
                             if(showRunningCircAvg)
                                 pR=plot(runningAvgSpeedBinCentersComp,phaseTimeCompRunningAvgValues,'r-','LineWidth',2);
                                 pR.Color=[1 0 0 runningAvgAlpha]
                             end
                             if(showLinCircFit)
                                 
                                 %validIdxesComp=phaseTimeCompRunningAvgValues<=maxCompressionVal;
                                 validIdxesComp=allUsedPhaseTimeCompressions<=maxCompressionVal;
                                 
                                %[m,b,R,p]=getLinearFit(runningAvgSpeedBinCentersComp(validIdxesComp),phaseTimeCompRunningAvgValues(validIdxesComp),currMinSpeed-0.2*currSpeedRange,currMaxSpeed+0.2*currSpeedRange)
                                 [m,b,R,p]=getLinearFit(allUsedRunSpeedFieldsPerSec(validIdxesComp),allUsedPhaseTimeCompressions(validIdxesComp),currMinSpeed-0.2*currSpeedRange,currMaxSpeed+0.2*currSpeedRange)

                                
                                minY=0;
                                maxY=max(allUsedPhaseTimeCompressions(validIdxesComp));
                                
                                try
                                 ylim([minY maxY])
                                end
                             end
                                 %{
                             legend(pR,'running avg')
                             legend boxoff
                                 %}
                             
                             if(currMaxDispSpeed>currMinDispSpeed)
                                xlim([currMinDispSpeed currMaxDispSpeed])
                          end
                             hold off
                         end 





                    uberTitle(removeUnderscores(sprintf('%s, %s',currSessName,seqStrs{ti})))
                    setFigFontTo(16)

                    %if(enoughPtsInEachCategory && highSpeedOccCount>=0 && lowSpeedOccCount>=0)
                    if(highSpeedOccCount>=0 && lowSpeedOccCount>=0)
                        try
                            maxFig
                            saveas(gcf,fullfile(exampleThetaSeqVsSpeedDir,sprintf('thetaSequencesVsSpeed_%s_%s.png',currSessName,seqStrs{ti})))
                        catch
                            maxFigMukkaalWidth
                            saveas(gcf,fullfile(exampleThetaSeqVsSpeedDir,sprintf('thetaSequencesVsSpeed_%s_%s.png',currSessName,seqStrs{ti})))
                        end
                    end
                    close all
                end
                
               
                if(showPlots && showSpeedSortedRaster)
                    
                    %subplot(1,2,2)
                    %subplot(2,2,[2 4])

                    for oi=1:numOccurences

                        currOccAllPhaseSeq=allPhaseSeqPerSpeedSortedOcc{oi};
                        currOccMeanPhaseSeq=meanPhaseSeqPerSpeedSortedOcc{oi};

                        %rasterAlpha=0.3;
                        rasterAlpha=0.4;
                        for ofi=1:length(currOccAllPhaseSeq)
                            currOccFieldPhases=currOccAllPhaseSeq{ofi};
                            pRaster=plotRasterStyle(currOccFieldPhases,oi,NaN,NaN,[distColors(fieldIdxToCenterRank(ofi),:) rasterAlpha]);
                            hold on
                            %show repeat cycle to emphasize circular nature
                            plotRasterStyle(currOccFieldPhases+360,oi,NaN,NaN,[distColors(fieldIdxToCenterRank(ofi),:) rasterAlpha]);
                            
                            
                            %pRaster.Color(4)=rasterAlpha;
                            %plot(currOccMeanPhaseSeq(ofi),oi,'.','Color',distColors(fieldIdxToCenterRank(ofi),:),'MarkerSize',round(20*20/numOccurences))
                            %plot(currOccMeanPhaseSeq(ofi),oi,'.','MarkerFaceColor',distColors(fieldIdxToCenterRank(ofi),:),'MarkerEdgeColor','k','MarkerSize',round(20*20/numOccurences))
                            %scatter(currOccMeanPhaseSeq(ofi),oi,round(50*20/numOccurences),'MarkerFaceColor',distColors(fieldIdxToCenterRank(ofi),:),'MarkerEdgeColor','k','LineWidth',1)
                            %scatter(currOccMeanPhaseSeq(ofi),oi,round(200*20/numOccurences),'MarkerFaceColor',distColors(fieldIdxToCenterRank(ofi),:),'MarkerEdgeColor','k','LineWidth',1)
                             scatter(currOccMeanPhaseSeq(ofi),oi,round(200*20/numOccurences),'MarkerFaceColor',distColors(fieldIdxToCenterRank(ofi),:),'MarkerEdgeColor',distColors(fieldIdxToCenterRank(ofi),:),'LineWidth',1)

                             
                             %show repeat cycle to emphasize circular nature
                             
                             scatter(currOccMeanPhaseSeq(ofi)+360,oi,round(200*20/numOccurences),'MarkerFaceColor',distColors(fieldIdxToCenterRank(ofi),:),'MarkerEdgeColor',distColors(fieldIdxToCenterRank(ofi),:),'LineWidth',1)

                             
                            %scatter(x,y,sz,'MarkerEdgeColor',[0 .5 .5],...
                            %  'MarkerFaceColor',[0 .7 .7],...
                            %  'LineWidth',1.5)
                            %plot(currOccMeanPhaseSeq(ofi),oi,'ko','MarkerSize',round(100/numOccurences))

                        end
                    end
                    
                     if(lowSpeedMaxCutoff>=min(sortedSpeedPerOcc))
                        [~,lowSpeedMaxCutoffIdx]=min(abs(sortedSpeedPerOcc-lowSpeedMaxCutoff));
                        lowSpeedMaxCutoffIdx=round(lowSpeedMaxCutoffIdx)+0.5;
                    else
                        lowSpeedMaxCutoffIdx=0;
                    end
                    
                    if(highSpeedMinCutoff<=max(sortedSpeedPerOcc))
                        [~,highSpeedMinCutoffIdx]=min(abs(sortedSpeedPerOcc-highSpeedMinCutoff));
                        highSpeedMinCutoffIdx=round(highSpeedMinCutoffIdx)-0.5;
                    else
                        highSpeedMinCutoffIdx=length(sortedSpeedPerOcc)+0.5;
                    end
                    
                    pL=hline(lowSpeedMaxCutoffIdx,'k--',3);
                    pH=hline(highSpeedMinCutoffIdx,'k-',3);
                    
                    %legend([pL pH],{sprintf('%.2f fields/sec',lowSpeedMaxCutoff),sprintf('%.2f fields/sec',highSpeedMinCutoff)});%,'Location','eastoutside')
                    %legend boxoff
                    box off
                    %xlim([0 360])
                    xlim([0 720])
                    ylim([0 numOccurences+1])
                    
                      xlabel('Theta phase (deg)')
                    ylabel('Speed-sorted traversal no.')
                    
                        %colormap(jet(numFieldsInSeq))
                        colormap(distColors)
                   cb=colorbar;
                   cb.Ticks=[]
                    ylabel(cb,'Behavioral field order')
                    caxis([1 numFieldsInSeq])
                    
                    
                
                
                uberTitle(removeUnderscores(sprintf('%s, %s',currSessName,seqStrs{ti})))
                    setFigFontTo(18)

                    %if(enoughPtsInEachCategory && highSpeedOccCount>=0 && lowSpeedOccCount>=0)
                    if(highSpeedOccCount>=0 && lowSpeedOccCount>=0)
                        
                            maxFig
                            saveas(gcf,fullfile(exampleThetaSeqVsSpeedDir,sprintf('thetaSequenceRastersVsSpeed_%s_%s.png',currSessName,seqStrs{ti})))

                    end
                end
                continue
                
                

            end

        
            
            if(false && showPlots && enoughPtsInEachCategory)
                subplot(1,2,1)
                %travSpeedPerOcc=cell2num(travAvgSpeedPerOcc);
                plot(travSpeedPerOcc,avgPhaseDiffPerOcc,'ko')

                hold on
                plot(travSpeedPerOcc(highSpeedIdxes),avgPhaseDiffPerOcc(highSpeedIdxes),'r.','MarkerSize',10)
                plot(travSpeedPerOcc(lowSpeedIdxes),avgPhaseDiffPerOcc(lowSpeedIdxes),'b.','MarkerSize',10)

                title('Avg phase diff in theta sequence vs traversal speed')

                 subplot(1,2,2)
                
                plot(travSpeedPerOcc,avgPhasePerOcc,'ko')
                hold on
                hold on
                plot(travSpeedPerOcc(highSpeedIdxes),avgPhasePerOcc(highSpeedIdxes),'r.','MarkerSize',10)
                plot(travSpeedPerOcc(lowSpeedIdxes),avgPhasePerOcc(lowSpeedIdxes),'b.','MarkerSize',10)
                
                title('Avg phase in theta sequence vs traversal speed')
                drawnow
                circFigH=fH;
            end
            
            if(false && showPlots && enoughPtsInEachCategory)
                subplot(1,2,1)
                box off
                 xlabel('Avg speed in traversal (m/s)')
                ylabel('Avg phase diff. in theta seq. (deg)')
                axis square
                
                %figure; 
                %close all
                 numImages=length(currSeqData.currThetaSeqSummaryImgPaths);
                 highSpeedImgCount=1;
                 lowSpeedImgCount=1;
                for ai=1:numImages
                    
                    if(highSpeedIdxes(ai))
                            figure(11)
                        subplot(floor(sqrt(sum(highSpeedIdxes))),floor(sqrt(sum(highSpeedIdxes)))+1,highSpeedImgCount)

                       imshow(currSeqData.currThetaSeqSummaryImgPaths{ai})
                       
                       highSpeedImgCount=highSpeedImgCount+1;
                   title(sprintf('Cycle %d',currSeqData.cycleNumPerOcc{ai}))
                    end
                    
                    if(lowSpeedIdxes(ai))
                            figure(12)
                        subplot(floor(sqrt(sum(lowSpeedIdxes))),floor(sqrt(sum(lowSpeedIdxes)))+1,lowSpeedImgCount)

                       imshow(currSeqData.currThetaSeqSummaryImgPaths{ai})
                      
                       
                       lowSpeedImgCount=lowSpeedImgCount+1;
                       title(sprintf('Cycle %d',currSeqData.cycleNumPerOcc{ai}))
                   
                    end
                    
                     
                    %autoArrangeFigures()
                    
                    
                end
                
                close all
            end
            
            
            
            circFigH=[]; %no circ corr lines
           
           %{ 
            [currSeqPhaseDiffSpeedRho,currSeqPhaseDiffSpeedPval,currSeqPhaseDiffSpeedSlopeDegPerMsec, ~]...
                =getCircCorrCoeff(travSpeedPerOcc(computeCorrIdxes),avgPhaseDiffPerOcc(computeCorrIdxes),circFigH);
            
            %}
            
            if(false && showPlots && enoughPtsInEachCategory)
                
                if(useCircMeasureSpread)
                    ylim([0 1])
                else
                    ylim([-180 180])
                end
            end
            
            
         
            if(false && showPlots && enoughPtsInEachCategory)
                subplot(1,2,2)
                box off
                xlabel('Avg speed in traversal (m/s)')
                ylabel('Avg phase of theta seq. (deg)')
                axis square
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %plot raster examples similar to Fig 1B
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if(numFieldsInSeq>=4)
                    %{
                    figure;
                    speedRank
                    
                    for oi=1:numOccurences
                        peakPhaseSeqPerOcc{oi}
                    end
                    %}
                    
                end
            end
            %{
            [currSeqMeanPhaseSpeedRho,currSeqMeanPhaseSpeedPval,currSeqMeanPhaseSpeedSlopeDegPerMsec, ~]...
                =getCircCorrCoeff(travSpeedPerOcc(computeCorrIdxes),avgPhasePerOcc(computeCorrIdxes),circFigH);
            %}
            
            
            
            
              if(false && showPlots && enoughPtsInEachCategory)
                ylim([0 360])
              end
            
              if(false && showPlots && enoughPtsInEachCategory)
                  uberTitle(removeUnderscores(sprintf('%s, %s',currSessName,seqStrs{ti})))
          
                  maxFigHalfHalfWidth
                  setFigFontTo(16)
                  saveas(gcf,fullfile(offsetAndDiffVsSpeedDir,sprintf('indivThetaSeqOffsetAndDiffVsSpeed_%s_Seq%s.png',currSessName,seqStrs{ti})))
                   cla
                  %close all
               
                  subplot(1,2,1)
                hold off
                subplot(1,2,2)
                hold off
               
              end
            
           
                
              %{
             if(abs(currSeqPhaseDiffSpeedSlopeDegPerMsec)<= maxAllowableSlopeMagnitude ...
                     && abs(currSeqMeanPhaseSpeedSlopeDegPerMsec)<=maxAllowableSlopeMagnitude)
                 
              %}
              %{
                 masterPhaseDiffSpeedRhos=[masterPhaseDiffSpeedRhos; currSeqPhaseDiffSpeedRho(:)];
                 masterMeanPhaseSpeedRho=[masterMeanPhaseSpeedRho; currSeqMeanPhaseSpeedRho(:)];
          
                 masterPhaseDiffSpeedSlope=[masterPhaseDiffSpeedSlope; currSeqPhaseDiffSpeedSlopeDegPerMsec(:)];
                 masterMeanPhaseSpeedSlope=[masterMeanPhaseSpeedSlope; currSeqMeanPhaseSpeedSlopeDegPerMsec(:)];
             %}

                 if(zScoreSpeed)
                   masterTravSpeed=[masterTravSpeed; zscoreLFP(travSpeedPerOcc(:))];
                 else
                    
                     masterTravSpeed=[masterTravSpeed; travSpeedPerOcc(:)];
                 end
                
                masterPhaseDiff=[masterPhaseDiff; avgPhaseDiffPerOcc(:)];
                masterMeanPhase=[masterMeanPhase; avgPhasePerOcc(:)];
            
             %end
             
             if(sum(highSpeedIdxes)<minNumPtsInSpeedCategory || sum(lowSpeedIdxes)<minNumPtsInSpeedCategory)
                %continue
             end
            
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %get high speed and low speed values
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            meanDistFracHighSpeed=nanmean(highSpeedDistFracs);
             meanDistFracLowSpeed=nanmean(lowSpeedDistFracs);
            
            
            meanPhaseDiffHighSpeed = circMeanDeg(avgPhaseDiffPerOcc(highSpeedIdxes));
             meanPhaseDiffLowSpeed = circMeanDeg(avgPhaseDiffPerOcc(lowSpeedIdxes));
             
              meanPhaseOffsetHighSpeed = circMeanDeg(avgPhasePerOcc(highSpeedIdxes));
               meanPhaseOffsetLowSpeed = circMeanDeg(avgPhasePerOcc(lowSpeedIdxes));
            
               if(nonCircStats)
                   masterPhaseDiffHighMinusLowSpeed=[masterPhaseDiffHighMinusLowSpeed;(meanPhaseDiffHighSpeed- meanPhaseDiffLowSpeed) ];
                 masterMeanPhaseHighMinusLowSpeed=[masterMeanPhaseHighMinusLowSpeed; (meanPhaseOffsetHighSpeed-meanPhaseOffsetLowSpeed) ];

               else
                 masterPhaseDiffHighMinusLowSpeed=[masterPhaseDiffHighMinusLowSpeed; angdiffDeg([meanPhaseDiffLowSpeed meanPhaseDiffHighSpeed])];
                 masterMeanPhaseHighMinusLowSpeed=[masterMeanPhaseHighMinusLowSpeed; angdiffDeg([meanPhaseOffsetLowSpeed meanPhaseOffsetHighSpeed])];
               end
            
            
            masterMeanDistFracHighMinusLowSpeed=[masterMeanDistFracHighMinusLowSpeed; (meanDistFracHighSpeed-meanDistFracLowSpeed)];
            %{
            if(showPlots)
                 figure; 
                 numImages=length(currSeqData.currThetaSeqSummaryImgPaths);
                for ai=1:numImages
                    subplot(floor(sqrt(numImages)),floor(sqrt(numImages))+1,ai)
                   imshow(currSeqData.currThetaSeqSummaryImgPaths{ai})
                end
                maxFig
                close all
            end
            %}
        end
        
        
    end
    
    %currThetaSeqData.thetaSequenceInfoPerSequence=thetaSequenceInfoPerSequence; %update sequence info
    
    save(currSeqDataFilePath,'thetaSequenceInfoPerSequence','-append')
    
end

                     
save('thetaSequencePropVsSpeedProcessedData.mat',...
    'meanDistPerFieldSequence', 'distSepPerFieldSequence','timeSepPerFieldSequence','runSpeedsPerFieldSequence','phaseSepPerFieldSequence','meanPhasePerFieldSequence','phaseCompPerFieldSequence','fieldSeqCount',...
    'allSessionRunSpeed','allSessionMeanDistFracs','allSessionAvgDistSepFracs','allSessionAvgTimeSepFracs','allSessionPhaseSep','allSessionMeanPhase','allSessionPhaseComp')
%%
%{

fH=figure
phaseDiffEdges=linspace(-180, 180, 71);
circMeasureSpreadEdges=linspace(0,0.2,71);

if(zScoreSpeed)
    %travSpeedEdges=linspace(-3,3,21);
    travSpeedEdges=linspace(minZspeed,maxZspeed,21);
 
else
    travSpeedEdges=linspace(minSpeed,maxSpeed,21);
end
phaseEdges=linspace(0,360,71);

withSmooth=1;

subplot(1,2,1)
%{
if(useCircMeasureSpread)
    masterPhaseDiff(masterPhaseDiff<0.0001)=NaN;
    getJointDistrGivenX(masterTravSpeed,masterPhaseDiff,travSpeedEdges,circMeasureSpreadEdges,fH,withSmooth)
else
    getJointDistrGivenX(masterTravSpeed,masterPhaseDiff,travSpeedEdges,phaseDiffEdges,fH,withSmooth)

end
%}
getJointDistrGivenX(masterTravSpeed,masterPhaseDiff,travSpeedEdges,phaseDiffEdges,fH,withSmooth)
%getJointDistrGivenXpeakNorm(masterTravSpeed,masterPhaseDiff,travSpeedEdges,phaseDiffEdges,fH,withSmooth)

xlabel('Traversal speed (m/s)')
 if(speedInFieldsPerSec) 
     xlabel('Traversal speed (fields/sec)')
 end
ylabel('Avg phase diff in theta sequence (deg)')
axis square
title(sprintf('Speed vs phase separation, n=%d theta sequences',length(masterMeanPhase)))


subplot(1,2,2)
getJointDistrGivenX(masterTravSpeed,masterMeanPhase,travSpeedEdges,phaseEdges,fH,withSmooth)
%getJointDistrGivenXpeakNorm(masterTravSpeed,masterMeanPhase,travSpeedEdges,phaseEdges,fH,withSmooth)

axis square

xlabel('Traversal speed (m/s)')
 if(speedInFieldsPerSec) 
     xlabel('Traversal speed (fields/sec)')
 end
ylabel('Avg phase in theta sequence (deg)')
title(sprintf('Speed vs phase offset, n=%d theta sequences',length(masterMeanPhase)))
setFigFontTo(16)
maxFig
saveas(gcf,'SpeedVsThetaSequencePhaseDiffAndOffset.png')
%%
%getBinnedCircAverages
if(zScoreSpeed)
  %travSpeedAvgEdges=linspace(-3,3,16);  
  travSpeedAvgEdges=linspace(minZspeed,maxZspeed,16);  
else
travSpeedAvgEdges=linspace(minSpeed,maxSpeed,16);
end
%travSpeedAvgEdges=linspace(minSpeed,maxSpeed,15);
%[speedBinCenters,phaseAvgValues] = getBinnedCircAverages(masterTravSpeed,masterMeanPhase,travSpeedAvgEdges);
%[speedBinCenters,phaseDiffValues] = getBinnedCircAverages(masterTravSpeed,masterPhaseDiff,travSpeedAvgEdges);
[speedBinCenters,phaseAvgValues] = getBinnedCircAverages(masterTravSpeed,masterMeanPhase,travSpeedAvgEdges);

travSpeedAvgDiffEdges=linspace(minSpeed,maxSpeed,16);

[speedBinDiffCenters,phaseDiffValues] = getBinnedCircAverages(masterTravSpeed,masterPhaseDiff,travSpeedAvgDiffEdges);



figure;
subplot(1,2,1)
plot(speedBinCenters,phaseAvgValues,'k.','MarkerSize',30)
xlabel('Traversal speed (m/s)')

if(speedInFieldsPerSec)
    xlabel('Traversal speed (fields/sec)')
end
ylabel('Avg phase in theta sequence (deg)')

ylim([140 240])
if(zScoreSpeed)
    xlim([minZspeed maxZspeed])
else
    xlim([minSpeed maxSpeed])
end
axis square
box off

subplot(1,2,2)
plot(speedBinDiffCenters,phaseDiffValues,'k.','MarkerSize',30)
xlabel('Traversal speed (m/s)')
if(speedInFieldsPerSec)
    xlabel('Traversal speed (fields/sec)')
end
ylabel('Avg phase diff in theta sequence (deg)')
ylim([0 100])
if(zScoreSpeed)
    xlim([minZspeed maxZspeed])
else
    xlim([minSpeed maxSpeed])
end
axis square
box off

saveas(gcf,fullfile(offsetAndDiffVsSpeedDir,'SpeedVsThetaSequencePhaseDiffAndOffset_BinnedAvgs.png'))
%figure; histogram(masterPhaseDiffSpeedRhos,100)
%figure; histogram(masterMeanPhaseSpeedRho,100)

%%

if(isempty(masterPhaseDiffHighMinusLowSpeed))
    
end 

[pMeanDiff,h,statsMeanDiff] = signrank(masterPhaseDiffHighMinusLowSpeed);

[pOffset,h,statsOffset] = signrank(masterMeanPhaseHighMinusLowSpeed);

gaussKernelWidthDeg=10;
gaussKernelWidthDeg=15;
%gaussKernelWidthDeg=5;
figure; 
phaseDiffEdges=linspace(-180,180,101);
phaseDiffBins=edgesToBins(phaseDiffEdges);
subplot(1,2,2);
%histogram(masterPhaseDiffHighMinusLowSpeed,phaseDiffEdges)

[pDistPhaseDiff,~]=getCircKernelDistr(masterPhaseDiffHighMinusLowSpeed,phaseDiffEdges,gaussKernelWidthDeg);
plot(phaseDiffBins,pDistPhaseDiff,'k-','LineWidth',4)
xlim([-180 180])
ylim([0 max(pDistPhaseDiff)])
%ylim([0 0.027])
hold on
xlabel('high speed circ-minus low speed phase difference (deg)')
ylabel('Probability')
vline(0,'k--',3)
title({'Speed vs phase separation of individual field sequence representations',sprintf('Wilcoxon signed rank test, z=%.2f, p=%.5f',statsMeanDiff.zval,pMeanDiff)})

box off
subplot(1,2,1);
%histogram(masterMeanPhaseHighMinusLowSpeed,phaseDiffEdges)
[pDistPhaseOffset,~]=getCircKernelDistr(masterMeanPhaseHighMinusLowSpeed,phaseDiffEdges,gaussKernelWidthDeg);
%plot(phaseDiffBins,pDistPhaseOffset,'k-','LineWidth',4)
plot(phaseDiffBins,pDistPhaseOffset,'r-','LineWidth',4)
xlim([-180 180])

ylim([0 max(pDistPhaseOffset)])


xlabel('high speed circ-minus low speed mean phase (deg)')
ylabel('Probability')
hold on
vline(0,'k--',3)
vline(circMeanDegNoMod360(masterMeanPhaseHighMinusLowSpeed(:)),'r--',3)
title({'Speed vs phase offset of individual field sequence representations',sprintf('Wilcoxon signed rank test, z=%.2f, p=%.5f',statsOffset.zval,pOffset)})
setFigFontTo(18)
maxFigHalfWidth
box off
maxFigMukkaalWidth
saveas(gcf,fullfile(offsetAndDiffVsSpeedDir,'indThetaSeqPhaseOffsetVsPhaseDiffsStats.png'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%no significant difference in field fraction between high and low speed gps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; histogram(masterMeanDistFracHighMinusLowSpeed,21)

%}
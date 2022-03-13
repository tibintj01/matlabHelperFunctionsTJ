close all; clear all; clc
%dataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles';

processedDataDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData';
unitDataDir=fullfile(processedDataDir,'unitSpikeInfoStructs');
dataDir=unitDataDir;

%perCycleDataDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData/perCycleUnitSeqData';
perCycleDataDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData/perFieldPerCycleData';

earlyVsLateCyclePairImgDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData/earlyVsLateCyclePairImgDir';
touchDir(earlyVsLateCyclePairImgDir)

%cycleDataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/';
showPlots=1;
showPlots=0;
zScorePlot=1;
%zScorePlot=0;
usePhaseDiff=1;

compareSpeedAndTimeDiff=1;
compareSpeedAndTimeDiff=0;
minSpeed=0.05;

cycleHalfDuration=0.06; %60 msec
cycleDur=0.12;


plotVert=1;

minNumSamplesPerPair=5;
%minNumSamplesPerPair=10;

maxSpaceDiff=1;
maxSpaceDiff=0.5;
%maxSpaceDiff=0.3;
maxSpaceDiff=1;
spaceFracMax=0.3*0.3;
spaceFracMax=0.15;
spaceFracMax=0.5;
spaceFracMax=0.2;
spaceFracMax=0.5;
spaceFracMax=2;
spaceFracMax=0.2;
spaceFracMax=5;

spaceFracMax=0.1;
%spaceFracMax=10;
%spaceFracMax=0.2;
%spaceFracMax=0.2;
%spaceFracMax=0.05;

%spaceFracMax=0.2;
%spaceFracMax=0.15;
%spaceFracMax=0.15;
%spaceFracMin=0.05;
spaceFracMin=0;
%spaceFracMin=0.01;
%spaceFracMin=0;
%spaceFracMin=0.05;
%spaceFracMax=0.05;

%ABSOLUTE_THRESH_NUM_CYCLES_COFIRING=10
ABSOLUTE_THRESH_NUM_CYCLES_COFIRING=20
ABSOLUTE_THRESH_NUM_CYCLES_COFIRING=100
ABSOLUTE_THRESH_NUM_CYCLES_COFIRING=50
ABSOLUTE_THRESH_NUM_CYCLES_COFIRING=30 %enough to estimate mean diff...
ABSOLUTE_THRESH_NUM_CYCLES_COFIRING=0 %enough to estimate mean diff...

maxSEMDiff=0.005; %5 msec


maxSEMDiff=0.002; %2 msec
maxSEMDiff=0.003; %3 msec
maxSEMDiff=0.0025; %2.5 msec 
maxSEMDiffFrac=0.5; %at most 50 percent of mean - p=0.0029, 376 pairs
maxSEMDiffFrac=0.333333; %at most 1/3 of mean -

maxWidthRatio=100;
%maxSEMDiffFrac=0.2; %at most 1/3 of mean -
%{

minSpeedForThetaDiff=0.2;
maxSpeedForThetaDiff=0.6;
%}

%minSpeedForThetaDiff=0;
%minSpeedForThetaDiff=0.1;
minSpeedForThetaDiff=0.05;

maxSpeedForThetaDiff=10;

lastFrac=0.35;
lastFrac=0.33333;
%lastFrac=0.4;

%lastFrac=0.4;
%lastFrac=0.35;
lastFrac=0.33333;
lastFrac=0.5;
%lastFrac=0.33333;

%lastFrac=0.2
cofiringThresh=0.2;
cofiringThresh=0.1;
%cofiringThresh=0.08;

lateCycleTimeThresh=cycleDur*(1-lastFrac);
earlyCycleTimeThresh=cycleDur*lastFrac;
numBins=30;
numBins=50;

%restrictToCorrectOrderCycles=0;
restrictToCorrectOrderCycles=1;

restrictToNonContainedFields=1;
%endDiffThreshFrac=0.1;
%startDiffThreshFrac=0.1;
endDiffThreshFrac=100;
startDiffThreshFrac=100;

useDiffMod360=0;

%{
lateCycleTimeThresh=cycleDur*0.6;
earlyCycleTimeThresh=cycleDur*0.4;
%}

%{
lateCycleTimeThresh=cycleDur*0.5;
earlyCycleTimeThresh=cycleDur*0.5;
%}

sqlDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/sqlQueries';

sqlSessionFiles=getRegexFileNames(sqlDir,'ec*csv');

numSessions=length(sqlSessionFiles);
sessionNames=cell(numSessions,1);

for si=1:numSessions
    sessionNames{si}=getSubstrBtwnSubstrs(sqlSessionFiles{si},'','_sql');
end

%sessionNames={'Achilles_11012013'};
%sessionNames={'Gatsby_08282013'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%loop through each session and load per cycle spike time data across units
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
masterEarlyDiffs=[];
masterLateDiffs=[];
numEarlyPhasePairs=0;
numLatePhasePairs=0;

allOverlappingPairCenterDiffs=[];

startSi=1;
%startSi=78;
%startSi=80;
%startSi=3;
 allSpeedVsTimeDiffR=[];
 

for si=startSi:length(sessionNames) 
    currSessName=sessionNames{si};
    currSeqDataFilePath=fullfile(perCycleDataDir,sprintf('perFieldPerCycleTimingData%s.mat',currSessName));
    
     thisSessionAllEarlyAndLateDiffs=[];
    
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
     
     numCyclesCofiring=zeros(numFields,numFields);
       
    numEarlyPhasePairs=0;
    numLatePhasePairs=0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %for each direction, loop through each pair of units and each cycle 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %for di=1:2

        timeDifferencesPerCycle=NaN(numFields,numFields,numCycles);
        phaseDiffPerCycle=NaN(numFields,numFields,numCycles);
        currTimeInMetaFieldPerCycle=NaN(numFields,numFields,numCycles);
        numCyclesThreshPerRow=NaN(numFields,1);
        numCyclesMaxPerRow=NaN(numFields,1);
        
        alreadyProcessed=false(numFields,numFields);
        
        currFieldCentersPerCycle=currThetaSeqData.fieldCenterPerCyclePerField;
        currFieldStartTimesPerCycle=currThetaSeqData.fieldEntryTimePerCyclePerField;
        currFieldEndTimesPerCycle=currThetaSeqData.fieldExitTimePerCyclePerField;
        
        
        allFieldIDsThisSession=currThetaSeqData.allFieldIDsThisSession;
        
         timeDiffVsTimeH=figure;
        %loop through all pairs of units
        numPairPlots=0;
        for field1Num=1:numFields
            for field2Num=1:numFields %vs upper triangle matrix
            %for unit2Num=(unit1Num+1):numUnits %vs upper triangle matrix

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %get timing difference within cycle, to compare with speed,
                %spatial difference, and time in meta-field
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if(alreadyProcessed(field1Num,field2Num) || field1Num==field2Num)
                    continue
                end
                
                field1ID=allFieldIDsThisSession(field1Num);
                [ui1,di1,dfi1]=recoverFieldInfoFromID(field1ID);
                
                field2ID=allFieldIDsThisSession(field2Num);
                [ui2,di2,dfi2]=recoverFieldInfoFromID(field2ID);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %if these fields come from the same unit, skip
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if(ui1==ui2) 
                    continue
                end
              
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %if these fields come from the different directions, skip
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if(di1~=di2) 
                    continue
                else
                    di=di1;
                end
                
                if(di==1)
                    currDiStr='rightward';
                else
                    currDiStr='leftward';
                end
                    
                for ci=1:numCycles
                    field1TimingCurrCycle=currThetaSeqData.avgSpikeTimesPerCyclePerField(ci,field1Num);
                    field2TimingCurrCycle=currThetaSeqData.avgSpikeTimesPerCyclePerField(ci,field2Num);
                    
                    
                    field1PhaseCurrCycle=currThetaSeqData.peakSpikePhasePerCyclePerField(ci,field1Num);
                    field2PhaseCurrCycle=currThetaSeqData.peakSpikePhasePerCyclePerField(ci,field2Num);
                    
                    %timeFromUnit1Start=timeOfCycleTrough(ci)-currFieldStartTimesPerCycle(ci,unit1Num);
                    %timeFromUnit2Start=timeOfCycleTrough(ci)-currFieldStartTimesPerCycle(ci,unit2Num);
                    
                    
                    timeFromUnit1Start=timeOfCycleTrough(ci)-currFieldStartTimesPerCycle(ci,field1Num);
                    timeFromUnit2Start=timeOfCycleTrough(ci)-currFieldStartTimesPerCycle(ci,field2Num);
                    
                    %currTimeInMetaFieldPerCycle(unit1Num,unit2Num,ci)=max(timeFromUnit1Start,timeFromUnit2Start);
                    currTimeInMetaFieldPerCycle(field1Num,field2Num,ci)=timeFromUnit2Start-timeFromUnit1Start;
                    
                    phaseDiffPerCycle(field1Num,field2Num,ci)=angdiffDeg([field1PhaseCurrCycle field2PhaseCurrCycle]);

                    timeDifferencesPerCycle(field1Num,field2Num,ci)=field2TimingCurrCycle-field1TimingCurrCycle;
                    %{
                    if(unit1Num~=unit2Num && abs(timeDifferencesPerCycle(unit1Num,unit2Num,ci))<1e-10)
                        disp('')
                        
                    end
                    %}
                    
                    if(~isnan(timeDifferencesPerCycle(field1Num,field2Num,ci)))
                        %if(abs(timeDifferencesPerCycle(field1Num,field2Num,ci))<cycleDur*lastFrac)
                            numCyclesCofiring(field1Num,field2Num)=numCyclesCofiring(field1Num,field2Num)+1;
                        %end
                    end
                end
                %{
                 timeDiffVsTime=squeeze(timeDifferencesPerCycle(unit1Num,unit2Num,:));
                figure(timeDiffVsTimeH)
                    plot(timeDiffVsTime(~isnan(timeDiffVsTime))+numPairPlots/10,'LineWidth',3)
                    hold on
                    if(sum(~isnan(timeDiffVsTime))>0)
                        numPairPlots=numPairPlots+1;
                    end
                %}
                    
                alreadyProcessed(field1Num,field2Num)=true;
                alreadyProcessed(field2Num,field1Num)=true;
                
            end
            field1Num
            %numCyclesCofiring(unit1Num,unit1Num,di)=NaN;
            numCyclesMaxPerRow(field1Num)=numCyclesCofiring(field1Num,field1Num);
            %numCyclesThreshPerRow(unit1Num)=0.1*numCyclesCofiring(unit1Num,unit1Num,di);
            %numCyclesThreshPerRow(unit1Num)=0.3*numCyclesCofiring(unit1Num,unit1Num,di);
            %numCyclesThreshPerRow(unit1Num)=0.2*numCyclesCofiring(unit1Num,unit1Num,di);
            %numCyclesThreshPerRow(unit1Num)=0.15*numCyclesCofiring(unit1Num,unit1Num,di);
            numCyclesThreshPerRow(field1Num)=cofiringThresh*numCyclesCofiring(field1Num,field1Num);
        end
        
        earlyThetaTimeDiffs=[];
       lateThetaTimeDiffs=[];
       
       allPlottedLateDiffs=[];
       numLatePhasePairs=0;
       numEarlyPhasePairs=0;

       speedPerCycle=currThetaSeqData.speedPerCycle;
       timeOfCycleTrough=currThetaSeqData.timeOfCycleTrough;
       
        for field1Num=1:numFields
            %for field2Num=1:numFields %vs upper triangle matrix  
            for field2Num=field1Num:numFields %vs upper triangle matrix  
                
                if(field1Num==field2Num)
                    continue
                end
                
                field1ID=allFieldIDsThisSession(field1Num);
                [ui1,di1,dfi1]=recoverFieldInfoFromID(field1ID);
                
                field2ID=allFieldIDsThisSession(field2Num);
                [ui2,di2,dfi2]=recoverFieldInfoFromID(field2ID);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %if these fields come from the same unit, skip
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if(ui1==ui2) 
                    continue
                end
              
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %if these fields come from the different directions, skip
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if(di1~=di2) 
                    continue
                else
                    di=di1;
                end
                
                if(di==1)
                    currDiStr='rightward';
                else
                    currDiStr='leftward';
                end
                
                %not overlapping
                %if(numCyclesCofiring(unit1Num,unit2Num,di) < numCyclesThreshPerRow(unit1Num))
                if(numCyclesCofiring(field1Num,field2Num) < ABSOLUTE_THRESH_NUM_CYCLES_COFIRING)
                    continue
                end
                
                %center1=nanmedian(currFieldCentersPerCycle(:,field1Num));
                %center2=nanmedian(currFieldCentersPerCycle(:,field2Num));
                center1=currThetaSeqData.fieldCenterPerField(field1Num);
                center2=currThetaSeqData.fieldCenterPerField(field2Num);
                
                width1=currThetaSeqData.fieldWidthPerField(field1Num);
                width2=currThetaSeqData.fieldWidthPerField(field2Num);

                biggerWidth=max([width1 width2]);
                smallerWidth=min([width1 width2]);
                
                endDiffThresh=smallerWidth*endDiffThreshFrac;
                startDiffThresh=smallerWidth*startDiffThreshFrac;
                
                widthRatio=biggerWidth/smallerWidth;
                
                if(widthRatio>maxWidthRatio)
                    continue
                end
                
                if(di==1) %rightward travel
                    if(center2<center1)
                        continue
                    end
                    
                    currPairOverlap=(currThetaSeqData.fieldEndPerField(field1Num)-currThetaSeqData.fieldStartPerField(field2Num));
                    
                    currPairEndDiff=(currThetaSeqData.fieldEndPerField(field2Num)-currThetaSeqData.fieldEndPerField(field1Num)); 
                    
                    currPairStartDiff=(currThetaSeqData.fieldStartPerField(field2Num)-currThetaSeqData.fieldStartPerField(field1Num)); 
                    
                else
                    if(center2>center1)
                        continue
                    end
                    currPairOverlap=-(currThetaSeqData.fieldEndPerField(field1Num)-currThetaSeqData.fieldStartPerField(field2Num));
                    
                    currPairEndDiff=-(currThetaSeqData.fieldEndPerField(field2Num)-currThetaSeqData.fieldEndPerField(field1Num));
                    
                    currPairStartDiff=-(currThetaSeqData.fieldStartPerField(field2Num)-currThetaSeqData.fieldStartPerField(field1Num)); 
                end
                
                if(currPairOverlap<0)
                        continue
                    end
                
                if(restrictToNonContainedFields)
                     if(currPairEndDiff<0)% || currPairEndDiff>endDiffThresh )
                         continue
                     end
                     
                     if(currPairStartDiff<0)% || currPairStartDiff>startDiffThresh )
                         continue
                     end
                end
                
                 fieldWidthSum=width1+width2;
                
                if(isnan(center1) || isnan(center2))
                    continue
                end
                 %currSpaceDiffMag=abs(center2-center1);
                 
                 %only consider pairs where 2 follows 1 in direction of
                 %motion
                 if(strcmp(currDiStr,'rightward'))
                  currSpaceDiffMag=(center2-center1);
                 else
                     currSpaceDiffMag=-(center2-center1);
                 end
                 
                 
                 
                 %currField1PosFracsPerCycle=currThetaSeqData.avgSpikeFieldFracPosPerCyclePerField{field1Num};
                  %   currField2PosFracsPerCycle=currThetaSeqData.avgSpikeFieldFracPosPerCyclePerField{field2Num};
                     
                     
                
                 %allOverlappingPairCenterDiffs=[allOverlappingPairCenterDiffs;currSpaceDiffMag]; 
                 
                 %maxSpaceDiff=1 in this case
                 %if(currSpaceDiffMag>maxSpaceDiff*spaceFracMax || currSpaceDiffMag<maxSpaceDiff*spaceFracMin) %ensures same directionality (2 beyond 1)
                  biggerWidth=max([width1 width2])
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 %if center difference is greater than half of bigger width, can expect
                 %this signal to span early and late phases
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 if(currSpaceDiffMag>width1/2 || currSpaceDiffMag>width2/2 || currSpaceDiffMag<0)
                     
                 %if(currSpaceDiffMag>biggerWidth/2 || currSpaceDiffMag<0)
                 %smallerWidth=min([width1 width2]);
                  %if(currPairOverlap<smallerWidth/2 || currSpaceDiffMag<0) %if overlap amount is less than half field for the smaller field, skip (not enough overlap)
                  %if(currSpaceDiffMag>maxSpaceDiff*spaceFracMax  )   
                 %if(currSpaceDiffMag>maxSpaceDiff )

                     continue
                 end
                 
                 if(~isnan(currSpaceDiffMag))
                   allOverlappingPairCenterDiffs=[allOverlappingPairCenterDiffs;currSpaceDiffMag]; 
                 end
                 %{
                  if(currSpaceDiffMag<maxSpaceDiff*spaceFracMin)
                        continue
                  end
                 %}
                 
                   thisCellEarlyPhaseDiffs=[];
                    thisCellLatePhaseDiffs=[];
                    
                    thisCellPhaseDiffPerCycle=NaN(numCycles,1);
                    %{
                    figure(101)
                    close(figure(101))
                    %}
                    
                    
                    field1SpikeTimePerCycle=NaN(numCycles,1);
                    field2SpikeTimePerCycle=NaN(numCycles,1);
                    
                    %avgLeadingFieldTiming
                for ci=1:numCycles
                    
                    if(abs(speedPerCycle(ci))<minSpeedForThetaDiff || speedPerCycle(ci)>maxSpeedForThetaDiff)
                        continue
                    end
                    
                    field1TimingCurrCycle=currThetaSeqData.avgSpikeTimesPerCyclePerField(ci,field1Num);
                    field2TimingCurrCycle=currThetaSeqData.avgSpikeTimesPerCyclePerField(ci,field2Num);
                    
                    field1PhaseCurrCycle=currThetaSeqData.peakSpikePhasePerCyclePerField(ci,field1Num);
                    field2PhaseCurrCycle=currThetaSeqData.peakSpikePhasePerCyclePerField(ci,field2Num);
                    
                    field1FieldFracCurrCycle=currThetaSeqData.avgSpikeFieldFracPosPerCyclePerField(ci,field1Num);
                    field2FieldFracCurrCycle=currThetaSeqData.avgSpikeFieldFracPosPerCyclePerField(ci,field2Num);
                    
                    timeCurrCycle=currThetaSeqData.timeOfCycleTrough(ci);
                    
                    
                    
                    
                    %field1PosFracCurrCycle=
                    %field2PosFracCurrCycle=
                    
                    if(isnan(field1PhaseCurrCycle) || isnan(field2PhaseCurrCycle))
                        continue
                    end
                    
                    currCycleStartTime=currThetaSeqData.timeOfCycleTrough(ci);
                    
                    %unit1SpikeTimePerCycle(ci)=unit1TimingCurrCycle+currCycleStartTime;
                    %unit2SpikeTimePerCycle(ci)=unit2SpikeTimePerCycle+currCycleStartTime;
                    
                    %currPhaseDiffMag=abs(timeDifferencesPerCycle(field1Num,field2Num,ci));
                    if(usePhaseDiff)
                        currPhaseDiffMag=phaseDiffPerCycle(field1Num,field2Num,ci);
                    else
                        currPhaseDiffMag=(timeDifferencesPerCycle(field1Num,field2Num,ci));
                    end
                    
                    
                     normalizeBySpaceDiff=0;
                if(normalizeBySpaceDiff)
                    currPhaseDiffMag=currPhaseDiffMag/currSpaceDiffMag;
                end
                
                
               
                
                
                    
                    thisCellPhaseDiffPerCycle(ci)=currPhaseDiffMag;
                    

                    
                    %if(field1Num~=field2Num && field1TimingCurrCycle > lateCycleTimeThresh && field2TimingCurrCycle >lateCycleTimeThresh)   
                    %if(field1Num~=field2Num && field1TimingCurrCycle > lateCycleTimeThresh && field2TimingCurrCycle >lateCycleTimeThresh)   
                    %if(field1PhaseCurrCycle > 180 && field2PhaseCurrCycle >180)   
                    if(field1PhaseCurrCycle > 360*(1-lastFrac) && field2PhaseCurrCycle >360*(1-lastFrac) && field1FieldFracCurrCycle<(1-lastFrac) && field2FieldFracCurrCycle<(1-lastFrac)) 

                        if(restrictToCorrectOrderCycles)
                            if(field1PhaseCurrCycle<field2PhaseCurrCycle)
                                
                        %lateThetaTimeDiffs=[lateThetaTimeDiffs; currPhaseDiffMag/currSpaceDiffMag];
                        %lateThetaTimeDiffs=[lateThetaTimeDiffs; currPhaseDiffMag];
                                thisCellLatePhaseDiffs=[thisCellLatePhaseDiffs; currPhaseDiffMag];
                            end
                        
                        else
                            thisCellLatePhaseDiffs=[thisCellLatePhaseDiffs; currPhaseDiffMag];
                        end
                        
                        if(isnan(currPhaseDiffMag))
                            %continue
                        end
                        
                    end
                    
                    %if(field1Num~=field2Num &&  field1TimingCurrCycle < earlyCycleTimeThresh && field2TimingCurrCycle <earlyCycleTimeThresh )
                    %if(field1PhaseCurrCycle < 180 && field2PhaseCurrCycle <180)  
                        
                        if(field1PhaseCurrCycle < 360*(lastFrac) && field2PhaseCurrCycle <360*(lastFrac) && field1FieldFracCurrCycle>lastFrac && field2FieldFracCurrCycle>lastFrac) 
                        %earlyThetaTimeDiffs=[earlyThetaTimeDiffs; currPhaseDiffMag/currSpaceDiffMag];
                        %earlyThetaTimeDiffs=[earlyThetaTimeDiffs; currPhaseDiffMag];
                        
                        if(restrictToCorrectOrderCycles)
                            if(field1PhaseCurrCycle<field2PhaseCurrCycle)
                            
                                thisCellEarlyPhaseDiffs=[thisCellEarlyPhaseDiffs; currPhaseDiffMag];
                            end
                        else
                            thisCellEarlyPhaseDiffs=[thisCellEarlyPhaseDiffs; currPhaseDiffMag];
                        end
                        if(isnan(currPhaseDiffMag))
                            %continue
                        end
                       
                    end
                    %drawnow
                end
                
                
                
                %{
                figure(112)
                plotRasterStyle(eventTimes,rowNum,startTick,endTick,tickColor)
                hold on
                plotRasterStyle(unit2SpikeTimePerCycle,rowNum,startTick,endTick,tickColor)
                %}
                
                %unit1EntryTimePerCycle=squeeze(currThetaSeqData.unitEntryTimePerCycleAndField(:,unit1Num,di));
                %{
                figure(111)
                plot(timeOfCycleTrough-unit1EntryTimePerCycle,thisCellPhaseDiffPerCycle)
                hold on
                xlim([0 3])
                drawnow
                %}
                
                if(isempty(thisCellLatePhaseDiffs) || isempty(thisCellEarlyPhaseDiffs))
                    continue
                end
                
               
                
                %if(getSEMacrossRows(thisCellLatePhaseDiffs)<=maxSEMDiff && getSEMacrossRows(thisCellEarlyPhaseDiffs)<=maxSEMDiff)
                %if(getSEMacrossRows(thisCellLatePhaseDiffs)<=(maxSEMDiffFrac*nanmean(thisCellLatePhaseDiffs)) && getSEMacrossRows(thisCellEarlyPhaseDiffs)<=(maxSEMDiffFrac*nanmean(thisCellEarlyPhaseDiffs)))

                if(length(thisCellEarlyPhaseDiffs)>=minNumSamplesPerPair && length(thisCellLatePhaseDiffs)>=minNumSamplesPerPair)
                    if(usePhaseDiff)
                        if(~useDiffMod360)
                            earlyThetaTimeDiffs=[earlyThetaTimeDiffs; circMeanDegNoMod360(thisCellEarlyPhaseDiffs)];
                            lateThetaTimeDiffs=[lateThetaTimeDiffs; circMeanDegNoMod360(thisCellLatePhaseDiffs)];
                        else
                             earlyThetaTimeDiffs=[earlyThetaTimeDiffs; circMeanDeg(thisCellEarlyPhaseDiffs)];
                            lateThetaTimeDiffs=[lateThetaTimeDiffs; circMeanDeg(thisCellLatePhaseDiffs)];
                        end
                        %earlyThetaTimeDiffs=[earlyThetaTimeDiffs; nanmean(thisCellEarlyPhaseDiffs)];
                        %lateThetaTimeDiffs=[lateThetaTimeDiffs; nanmean(thisCellLatePhaseDiffs)];
                    else
                    
                        earlyThetaTimeDiffs=[earlyThetaTimeDiffs; nanmean(thisCellEarlyPhaseDiffs)];
                        lateThetaTimeDiffs=[lateThetaTimeDiffs; nanmean(thisCellLatePhaseDiffs)];
                    end
                    disp('')
                else
                    continue
                end
                
                
                
                
                     %earlyThetaTimeDiffs=[earlyThetaTimeDiffs; nanmedian(thisCellEarlyPhaseDiffs)];
                    %lateThetaTimeDiffs=[lateThetaTimeDiffs; nanmedian(thisCellLatePhaseDiffs)];
                %end
                
                %{
                length(earlyThetaTimeDiffs)
                length(lateThetaTimeDiffs)
                
                if(length(lateThetaTimeDiffs)>175)
                disp('')
                end
                %}
                
                
                if(isnan(nanmean(thisCellEarlyPhaseDiffs)) || isnan(nanmean(thisCellLatePhaseDiffs)))
                    continue
                end
                %currRasterFig=figure(131+di)
                
                
                 currRasterFig=figure(131)
                 
                
                
                if(usePhaseDiff)
                    
                    if(~useDiffMod360)
                        earlyDiff=circMeanDegNoMod360(thisCellEarlyPhaseDiffs);
                        lateDiff=circMeanDegNoMod360(thisCellLatePhaseDiffs);
                    else
                    
                        earlyDiff=circMeanDeg(thisCellEarlyPhaseDiffs);
                        lateDiff=circMeanDeg(thisCellLatePhaseDiffs);
                    end
                    
                else
                    earlyDiff=nanmean(thisCellEarlyPhaseDiffs)*1000;
                    lateDiff=nanmean(thisCellLatePhaseDiffs)*1000;
                end
                
                if(exist('allPlottedLateDiffsSortedSaved','var'))
                    
                    %plotIdx=find(abs(lateDiff-allPlottedLateDiffsSortedSaved)<0.001);
                     %plotIdx=find(abs(abs(earlyDiff-lateDiff)-allPlottedLateDiffsSortedSaved)<0.001);
                     plotIdx=find(abs((earlyDiff+lateDiff)/2-allPlottedLateDiffsSortedSaved)<0.0001);
                else
                    plotIdx=numLatePhasePairs+1;
                end

                %subplot(2,1,2)
                        %plotRasterStyle(0,length(lateThetaTimeDiffs),NaN,NaN,NaN,3)
                     
                         numLatePhasePairs=numLatePhasePairs+1;
                        if(~plotVert)
                         pLate=plotRasterStyle(lateDiff,numLatePhasePairs,NaN,NaN,[1 0 0],3);
                        else
                            pLate=plotRasterStyleVert(lateDiff,plotIdx,NaN,NaN,[1 0 0],3)
                        end 
                        hold on
                        xlim([0 cycleDur*lastFrac*0.6]*1000)
                       
                        ylim([0 numLatePhasePairs])
                        
                 %subplot(2,1,1)
                        %plotRasterStyle(0,length(earlyThetaTimeDiffs),NaN,NaN,NaN,1)
                        %hold on
                        %pEarly=plotRasterStyle(earlyDiff,numEarlyPhasePairs,NaN,NaN,NaN,3)
                            numEarlyPhasePairs=numEarlyPhasePairs+1;
                        if(plotVert)
                            pEarly=plotRasterStyleVert(earlyDiff,plotIdx,NaN,NaN,NaN,3)
                        else
                            pEarly=plotRasterStyle(earlyDiff,numEarlyPhasePairs,NaN,NaN,NaN,3)

                        end
                      
                    
                        
                         if(~plotVert)
                             xlim([0 cycleDur*lastFrac*0.6]*1000)
                            ylim([0 numEarlyPhasePairs+1])
                            xlabel('Avg spike time difference (msec)')
                            ylabel('Cell pair no.')
                         else
                             if(usePhaseDiff)
                                 ylim([-180 180])
                                 ylim([0 130])
                                  ylabel('Avg spike phase difference (deg)')
                             else
                                ylim([0 cycleDur*lastFrac*0.6]*1000)
                                ylabel('Avg spike time difference (msec)')
                             end
                            xlim([0 numEarlyPhasePairs+1])
                            
                            xlabel('Cell pair no.')
                         end
                        
                        if(earlyDiff > lateDiff)
                            lineColor='b';
                        else
                            lineColor='g';
                        end
                        
                        if(~plotVert)
                            plot([earlyDiff lateDiff],[numEarlyPhasePairs numEarlyPhasePairs]-1,lineColor,'LineWidth',2)
                        else
                            plot([plotIdx plotIdx],[earlyDiff lateDiff],lineColor,'LineWidth',2)
                        end
                        hold on
                        
                        %allPlottedLateDiffs=[allPlottedLateDiffs; abs(earlyDiff-lateDiff)];
                        allPlottedLateDiffs=[allPlottedLateDiffs; (earlyDiff+lateDiff)/2];

                %{
                figure(101); hEarly=histogram(thisCellEarlyPhaseDiffs,100)
                hold on
                hLate=histogram(thisCellLatePhaseDiffs,100)
                hEarly.FaceColor='k';
                hLate.FaceColor='r';
                drawnow
                pause 0.1
                %}
                        
               thisSessionAllEarlyAndLateDiffs=[thisSessionAllEarlyAndLateDiffs ; [earlyDiff lateDiff]];
    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %plot current field pair for which early vs late
                %differences are calculated
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                currField1Phases=currThetaSeqData.currSessMasterSpikePhases{field1Num};
                currField1Positions=currThetaSeqData.currSessMasterSpikePositions{field1Num};
                
                currField2Phases=currThetaSeqData.currSessMasterSpikePhases{field2Num};
                currField2Positions=currThetaSeqData.currSessMasterSpikePositions{field2Num};
                
                 if(isempty(currField1Phases) || isempty(currField2Phases))
                     continue
                 end
                     
                 figure;
                        
                
                plot(currField1Positions,currField1Phases,'b.')
                hold on
                plot(currField2Positions,currField2Phases,'r.')
                legend('field 1', 'field 2')
                
                title({sprintf('Early cycle avg phase diff: %d deg', round(earlyDiff)),sprintf('Late cycle avg phase diff: %d deg', round(lateDiff))})
                figure(currRasterFig)
                
                clef
                        autoArrangeFigures();
                        
            end %field 2 loop
        end %field 1 loop
        
        
        clef
        saveAllOpenFigures(sprintf('allPairFields_%s.pdf',currSessName))
        
        [allPlottedLateDiffsSorted,sortByLateDiffIdx]=sort(allPlottedLateDiffs);
        %save(sprintf('plottedLateDiffs%s_%d',currSessName,di),'allPlottedLateDiffsSorted','sortByLateDiffIdx')
        
        if(isempty(thisSessionAllEarlyAndLateDiffs))
            continue
        end 
        
       setFigFontTo(18)
      % try
      
      try
           legend([pEarly pLate],'early cycle phase','late cycle phase')
           legend boxoff
      end
           

           if(usePhaseDiff)
                if(~plotVert)
                xlim([-30 200])
           else
                %ylim([-30 200])
                %ylim([0 120])
                 if(useDiffMod360)
                      ylim([0 360])
                 else
                     %ylim([-180 180])
                     ylim([0 130])
                 end
                    
           end
           else
           if(~plotVert)
                xlim([-15 36.5])
           else
                ylim([-15 36.5])
           end

           end


          box off
             %saveas(gcf,fullfile(earlyVsLateCyclePairImgDir,sprintf('%s_%sLateVsEarlyCyclePhaseDiffsPairwise.png',currSessName,currDiStr)))
             saveas(gcf,fullfile(earlyVsLateCyclePairImgDir,sprintf('%sLateVsEarlyCyclePhaseDiffsPairwise.png',currSessName)))
              saveEPS(fullfile(earlyVsLateCyclePairImgDir,sprintf('%sLateVsEarlyCyclePhaseDiffsPairwise.eps',currSessName)))
              
              figure;
              
              numPairs=size(thisSessionAllEarlyAndLateDiffs,1);
              for p=1:numPairs
                plot([2],thisSessionAllEarlyAndLateDiffs(p,1),'k.','MarkerSize',45*.8)
                hold on
                plot([1],thisSessionAllEarlyAndLateDiffs(p,2),'r.','MarkerSize',45*.8)
                
                if(thisSessionAllEarlyAndLateDiffs(p,2)<thisSessionAllEarlyAndLateDiffs(p,1))
                    plot([2 1],thisSessionAllEarlyAndLateDiffs(p,:),'b-','LineWidth',3)
                else
                    plot([2 1],thisSessionAllEarlyAndLateDiffs(p,:),'g-','LineWidth',3)
                end
              
              end
             
              xlim([0.75 2.25])
              
              ylim([20 105])
              xlabel('Theta phase of firing')
              ylabel('Phase difference of cell pair (deg)')
              setFigFontTo(18)
              box off
             %save( 'thisSessionAllEarlyAndLateDiffs','-append')
              %thisSessionAllEarlyAndLateDiffs=[];

              try
                 close(currRasterFig)
              end
        % end
        %figure; histogram(allOverlappingPairCenterDiffs,50)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %plot time differences vs speed, pos, time for each unique pair
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
        speedColors=magma(100);
        speedColors=parula(100);
        
        %edges=linspace(0,0.04,30); 
        %edges=linspace(0,0.05,100);
         %edges=linspace(0,0.05,numBins);
         
         if(usePhaseDiff)
             edges=linspace(0,360*lastFrac*0.6,numBins);
         else
            edges=linspace(0,cycleDur*lastFrac*0.6,numBins);
         end
         %edges=linspace(0,cycleDur*lastFrac*0.55,numBins);
         
         nonNaNIdxes=~isnan(earlyThetaTimeDiffs) & ~isnan(lateThetaTimeDiffs);
         
         numEarlyPairs=sum(~isnan(earlyThetaTimeDiffs(nonNaNIdxes)));
         numLatePairs=sum(~isnan(lateThetaTimeDiffs(nonNaNIdxes)));
         
         %{
         if(numEarlyPairs<15)
             continue
         end
         %}
         
         if(~isempty(earlyThetaTimeDiffs(nonNaNIdxes)))
            [p,h,stats] = signrank(earlyThetaTimeDiffs(nonNaNIdxes),lateThetaTimeDiffs(nonNaNIdxes));
            p
            if(isfield(stats,'zval'))
                zval=stats.zval
                signedrank=stats.signedrank
            end
            %{
            Ach10_25_13 rightward
            p = 4.5637e-04
              zval: 3.5051
             signedrank: 529
            
            Ach10_25_13 leftward
            p =0.0074
          zval: 2.6795
    signedrank: 281
            
            Ach11 leftward
            p =5.5166e-15
            zval =7.8145
            signedrank = 4225
            
            %}
         end
         
         masterEarlyDiffs=[masterEarlyDiffs; earlyThetaTimeDiffs(:)];
         masterLateDiffs=[masterLateDiffs; lateThetaTimeDiffs(:)];
         
         
        figure;
        subplot(2,1,1)
        
        imshow(currThetaSeqData.allFieldsInSessionImgPath)
        subplot(2,1,2)
        
        
        Ne=histcounts((earlyThetaTimeDiffs),edges); 
        Nl=histcounts((lateThetaTimeDiffs),edges); 
        phaseDiffBins=edgesToBins(edges);
        %kernelPhaseAxis=linspace(0,120,1000);
        if(useDiffMod360)
                   kernelPhaseAxis=linspace(0,360,1000);
         else
             kernelPhaseAxis=linspace(-180,180,1000);
         end
        %
   
        kernelWidth=5;
         kernelDistrEarlyPhaseDiffObj = fitdist(earlyThetaTimeDiffs(:),'Kernel','BandWidth',kernelWidth); %gaussian kernel estimation, 5 degree
         kernelDistrLatePhaseDiffObj = fitdist(lateThetaTimeDiffs(:),'Kernel','BandWidth',kernelWidth); %gaussian kernel estimation, 5 degree

                
                kernelDistrEarlyPhaseDiff=pdf(kernelDistrEarlyPhaseDiffObj,kernelPhaseAxis);
                kernelDistrLatePhaseDiff=pdf(kernelDistrLatePhaseDiffObj,kernelPhaseAxis);

        
        if(~usePhaseDiff)
            plot(1000*edgesToBins(edges),smooth(Ne/sum(Ne),5),'k-','LineWidth',5);
            hold on; 
            Nl=histcounts((lateThetaTimeDiffs),edges); 
            plot(1000*edgesToBins(edges),smooth(Nl/sum(Nl),5),'r-','LineWidth',4);
             xlabel('Avg spike time difference in theta cycle (msec)')
        else
            
           %plot(phaseDiffBins,smooth(Ne/sum(Ne),5),'k-','LineWidth',5);
           plot(kernelPhaseAxis,kernelDistrEarlyPhaseDiff,'k-','LineWidth',5);
            hold on; 
             plot(kernelPhaseAxis,kernelDistrLatePhaseDiff,'r-','LineWidth',5);
           % plot(phaseDiffBins,smooth(Nl/sum(Nl),5),'r-','LineWidth',4);
            
            
             xlabel('Avg spike phase difference in theta cycle (deg)')
        end
       
        
        ylabel('Probability')
        legend(sprintf('first %d%% of cycle (n=%d pairs)',round(lastFrac*100),numEarlyPairs),sprintf('last %d%% of cycle (n=%d pairs)',round(lastFrac*100),numLatePairs))
        %title(removeUnderscores(sprintf('%s, overlapping field pairs between %d and %d cm apart across all theta cycles, %s',currSessName,round(spaceFracMin*100),round(spaceFracMax*100),currDiStr)))
        %title(removeUnderscores(sprintf('%s, overlapping field pairs across all theta cycles, %s',currSessName,currDiStr)))
        title(removeUnderscores(sprintf('%s, overlapping field pairs across all theta cycles',currSessName)))


        setFigFontTo(18)
        legend boxoff
        box off
        thisSessionAllEarlyAndLateDiffs
        %saveas(gcf,fullfile(earlyVsLateCyclePairImgDir,sprintf('%s_%sLateVsEarlyCyclePhaseDiffs.png',currSessName,currDiStr)))
        saveas(gcf,fullfile(earlyVsLateCyclePairImgDir,sprintf('%sLateVsEarlyCyclePhaseDiffs.png',currSessName)))
        saveEPS(fullfile(earlyVsLateCyclePairImgDir,sprintf('%sLateVsEarlyCyclePhaseDiffs.png',currSessName)))
         
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %just late vs early cycle
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if(si==length(sessionNames))
            
             if(usePhaseDiff)
                 
                 %edges=linspace(0,120,50+1);
                 %edges=linspace(-180,180,100+1);
                  %edges=linspace(0,180,50+1);
                  if(useDiffMod360)
                    edges=linspace(0,360,100+1);
                  else
                      edges=linspace(-180,180,100+1);
                      edges=linspace(0,130,50+1);
                  end
            figure;Ne=histcounts((masterEarlyDiffs),edges); 
            plot(edgesToBins(edges),smooth(Ne/sum(Ne),5),'k-','LineWidth',5);
            hold on; 
            Nl=histcounts((masterLateDiffs),edges); 
            plot(edgesToBins(edges),smooth(Nl/sum(Nl),5),'r-','LineWidth',5);
                 
            %{
            kernelPhaseAxis=linspace(0,120,1000);
        kernelWidth=5;
         kernelDistrEarlyPhaseDiffObj = fitdist(masterEarlyDiffs(:),'Kernel','BandWidth',kernelWidth); %gaussian kernel estimation, 5 degree
         kernelDistrLatePhaseDiffObj = fitdist(masterLateDiffs(:),'Kernel','BandWidth',kernelWidth); %gaussian kernel estimation, 5 degree

                
                kernelDistrEarlyPhaseDiff=pdf(kernelDistrEarlyPhaseDiffObj,kernelPhaseAxis);
                kernelDistrLatePhaseDiff=pdf(kernelDistrLatePhaseDiffObj,kernelPhaseAxis);

                     plot(kernelPhaseAxis,kernelDistrEarlyPhaseDiff,'k-','LineWidth',5);
            hold on; 
             plot(kernelPhaseAxis,kernelDistrLatePhaseDiff,'r-','LineWidth',5);
            %}
            
            nonNaNIdxesMaster=~isnan(masterEarlyDiffs) & ~isnan(masterLateDiffs);
            
            
             else
                 figure;Ne=histcounts((masterEarlyDiffs),edges); 
            plot(1000*edgesToBins(edges),smooth(Ne/sum(Ne),5),'k-','LineWidth',5);
            hold on; 
            Nl=histcounts((masterLateDiffs),edges); 
            plot(1000*edgesToBins(edges),smooth(Nl/sum(Nl),5),'r-','LineWidth',5);
                 
             end
            
            [p,h,stats] = signrank(masterEarlyDiffs(nonNaNIdxesMaster),masterLateDiffs(nonNaNIdxesMaster))
            zval=stats.zval;
            %{
            p =

               9.2703e-17


            h =

              logical

               1


            stats = 

              struct with fields:

                      zval: 8.3138
                signedrank: 51433
            
            nanmean(masterEarlyDiffs(nonNaNIdxesMaster))*1000 =

                   16.5455 msec

           nanmean(masterLateDiffs(nonNaNIdxesMaster))*1000 =

                   15.9739 msec
            %}
            if(usePhaseDiff)
              xlabel('Avg spike phase difference in theta cycle (deg)')

            else
              xlabel('Avg spike time difference in theta cycle (msec)')

            end
            ylabel('Probability')
            legend(sprintf('first %d%% of cycle (n=%d pairs)',round(lastFrac*100),length(masterEarlyDiffs(nonNaNIdxesMaster))),sprintf('last %d%% of cycle (n=%d pairs)',round(lastFrac*100),length(masterLateDiffs(nonNaNIdxesMaster))))
            %title(removeUnderscores(sprintf('All overlapping field pairs between %d and %d cm apart across all theta cycles',round(spaceFracMin*100),round(spaceFracMax*100))))
            title(removeUnderscores(sprintf('Overlapping field pairs across all theta cycles')))
            setFigFontTo(18)
            legend boxoff
            saveas(gcf,fullfile(earlyVsLateCyclePairImgDir,sprintf('AllLateVsEarlyCyclePhaseDiffs.png')))
            
            saveEPS(fullfile(earlyVsLateCyclePairImgDir,sprintf('AllLateVsEarlyCyclePhaseDiffs')))

            if(usePhaseDiff)
                %diffDiffEdges=linspace(-70,70,100+1);
                
                 if(useDiffMod360)
                    diffDiffEdges=linspace(0,360,70+1);
                 else
                     diffDiffEdges=linspace(-180,180,70+1);
                 end
                 
                  diffDiffEdges=linspace(-70,70,40+1);
                %diffDiffEdges=linspace(-60,60,50+1);
                 figure;Ne=histcounts((masterEarlyDiffs-masterLateDiffs),diffDiffEdges); 
            plot(edgesToBins(diffDiffEdges),smooth(Ne/sum(Ne),4),'k-','LineWidth',5);
                        xlabel('Early - Late, avg spike phase difference in theta cycles (deg)')

            else
                 diffDiffEdges=linspace(-0.02,0.02,100+1);
                 figure;Ne=histcounts((masterEarlyDiffs-masterLateDiffs),diffDiffEdges); 
            plot(1000*edgesToBins(diffDiffEdges),smooth(Ne/sum(Ne),5),'k-','LineWidth',5);
                        xlabel('Early - Late, avg spike time difference in theta cycles (msec)')

            end
            
            hold on
            plot([0 0],ylim,'k--','LineWidth',3)

            ylabel('Probability')
             legend(sprintf('n=%d pairs',length(masterEarlyDiffs(nonNaNIdxesMaster))))

            %legend(sprintf('first %d%% of cycle (n=%d pairs)',round(lastFrac*100),numEarlyPairs),sprintf('last %d%% of cycle (n=%d pairs)',round(lastFrac*100),numLatePairs))
            %title(removeUnderscores(sprintf('%s, overlapping field pairs < %d cm apart across all theta cycles, %s',currSessName,round(spaceFracMax*100),currDiStr)))
            %title({removeUnderscores(sprintf('Overlapping field pairs across all theta cycles')),'Wilcoxon signed rank test: zval=8.3138, p=9.2703e-17'})
            title({removeUnderscores(sprintf('Overlapping field pairs across all theta cycles')),sprintf('Wilcoxon signed rank test: zval=%.3f, p=%.6f',zval,p)})

            %xlim([-70 70])
            setFigFontTo(18)
            legend boxoff
            box off
            saveas(gcf,fullfile(earlyVsLateCyclePairImgDir,sprintf('AllLateVsEarlyCyclePhaseDiffDiffsDistr.png')))
            saveEPS(fullfile(earlyVsLateCyclePairImgDir,sprintf('AllLateVsEarlyCyclePhaseDiffDiffsDistr')))
        end
        
       
        %continue
        
        %save(fullfile(processedDataDir,'allUnit
        
      if(compareSpeedAndTimeDiff)
        figure
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for fi=1:numFields
            for fii=(fi+1):numFields
                currPairTimeDiffsPerCycle=squeeze(timeDifferencesPerCycle(fi,fii,:));
                
                currFieldPairPhaseDiffPerCycle=phaseDiffPerCycle(fi,fii,:);
                
                
                %currPairTimeInMetaFieldPerCycle=squeeze(currTimeInMetaFieldPerCycle(ui,uii,:));
                
                center1PerCycle=currFieldCentersPerCycle(:,fi);
                center2PerCycle=currFieldCentersPerCycle(:,fii);
                
                if(numCyclesCofiring(fi,fii)<numCyclesThreshPerRow(fi) || numCyclesCofiring(fi,fii)<numCyclesThreshPerRow(fii))
                    continue
                end
                
                currPairSpaceDiffPerCycle=center2PerCycle-center1PerCycle;
                
                if(isnan(max(currPairTimeDiffsPerCycle)))
                    continue
                end
                
                validIdxes=~isnan(currPairTimeDiffsPerCycle) & abs(speedPerCycle)>minSpeed;
                currSpeedPerCycle=speedPerCycle(validIdxes);
                currPairTimeDiffsPerCycle=currPairTimeDiffsPerCycle(validIdxes);
                currPairSpaceDiffPerCycle=currPairSpaceDiffPerCycle(validIdxes);
                %currPairTimeInMetaFieldPerCycle=currPairTimeInMetaFieldPerCycle(validIdxes);
                
                thisFieldMeanSpeed=nanmean((currSpeedPerCycle));
                speedColorIdx=ceil(thisFieldMeanSpeed*100);
                currFieldPairMeanTimeDiff=nanmean(currPairTimeDiffsPerCycle);
                currFieldPairMeanSpaceDiff=nanmean(currPairSpaceDiffPerCycle);
                
                %subplot(2,1,1)
                if(zScorePlot)
                    %subplot(2,2,1)
                      subplot(2,1,1)
                    plot(zscoreLFP(currSpeedPerCycle),zscoreLFP(currPairTimeDiffsPerCycle),'k.')
                    xlabel('Running speed (Z)')
                    ylabel('Time difference within theta cycle (Z)')
                    xlim([-3 3])
                     ylim([-3 3])
                     axis square
                      hold on
                      
                     nonNanIdxDiffs=~isnan(currSpeedPerCycle) & ~isnan(currPairTimeDiffsPerCycle);
                     
                     if(length(currSpeedPerCycle(nonNanIdxDiffs))>1)
                         rMat=corrcoef(currSpeedPerCycle(nonNanIdxDiffs),currPairTimeDiffsPerCycle(nonNanIdxDiffs));
                         speedVsTimeDiffR=rMat(1,2);

                         allSpeedVsTimeDiffR=[allSpeedVsTimeDiffR; speedVsTimeDiffR];
                     end
                     %{
                     subplot(2,2,2)
                     histogram(allSpeedVsTimeDiffR,30)
                     
                     %}
                else
                    subplot(2,1,1)
                    plot(nanmean(currSpeedPerCycle),nanmean(currPairTimeDiffsPerCycle),'k.','MarkerSize',20)
                    xlabel('Avg running speed (m/s)')
                    ylabel('Avg time difference within theta cycle (sec)')
                     xlim([0.2 1])
                     if(contains(currSessName,'11012013'))
                        ylim([-0.04 0.03])
                     else
                          ylim([-0.06 0.06])
                     end
                      hold on
                 %ylim([-Inf Inf])     
                end
               
                
                subplot(2,1,2)
                 %plot(currPairSpaceDiffPerCycle,currPairTimeDiffsPerCycle,'k.')
                 %{
                 plot(thisUnitMeanSpaceDiff,thisUnitMeanTimeDiff,'.','MarkerSize',20,'Color',speedColors(speedColorIdx,:))
                 colorbar
                 caxis([0 1])
                 %}
                %scatter(thisUnitMeanSpaceDiff,thisUnitMeanTimeDiff,20,thisUnitMeanSpeed
                 
                 plot(currFieldPairMeanSpaceDiff,currFieldPairMeanTimeDiff,'k.','MarkerSize',20)
                    xlabel('Spatial difference between field centers (m)')
                    ylabel('Avg time difference within theta cycle (sec)')
                    
                hold on
                xlim([-0.5 0.5])
                     if(contains(currSessName,'11012013'))
                        ylim([-0.05 0.05])
                     else
                          ylim([-0.06 0.06])
                     end
                 %ylim([-Inf Inf])
                 
                 %{
                 subplot(3,1,3)
                  plot(nanmean(currPairTimeInMetaFieldPerCycle),nanmean(currPairTimeDiffsPerCycle),'k.','MarkerSize',20)
                  

                  xlabel('Time since field entry (sec)')
                    ylabel('Time difference within theta cycle (sec)')
                   xlim([-1 1])
                       %ylim([-0.06 0.06])
                     %axis square
                     hold on
                 %}
                 
               end
            end
        
        uberTitle(removeUnderscores(sprintf('%s, all overlapping field pairs across all theta cycles, %s',currSessName,currDiStr)))
                %maxFigFourthWidthFullHeight
                maxFigHalfWidth
                setFigFontTo(18)
                saveas(gcf,fullfile(earlyVsLateCyclePairImgDir,sprintf('%s_%sSpeedInvariantSpatialSeqCoding.png',currSessName,currDiStr)))
                
               
                %speedDiffRedges=linspace(-0.5,0.5,50);
                speedDiffRedges=linspace(-1,1,100);
                speedDiffBinCenters=edgesToBins(speedDiffRedges);
                %pSpeedThetaDiff=getProbDist(allSpeedVsTimeDiffR,speedDiffRedges,1,0);
                 pSpeedThetaDiff=getProbDist(allSpeedVsTimeDiffR,speedDiffRedges,0,0);
                
                
                 figure
                plot(speedDiffBinCenters,pSpeedThetaDiff,'k-')
                xlabel('Speed-phase difference correlation coefficient')
                legend(sprintf('n=%d field pairs (all sessions)',length(allSpeedVsTimeDiffR)))
                axis tight
        
                hold on; plot([0 0],ylim,'k--')
                ylabel('Probability')
                setFigFontTo(18)
                xlim([-1 1])
                        disp('')
                        box off
                        legend boxoff
                        maxFigHalfHalfWidth
             end
                
                %{
                nanmean(allSpeedVsTimeDiffR) = -0.0083
                getSEMacrossRows(allSpeedVsTimeDiffR) = 0.0085
                
                Signed rank test different from 0:
                p =0.3567
                zval: -0.9217
                signedrank: 21189
                
                %}
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %plot pairwise time differences heat map
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %{
        avgTimeDiffPerCycle=nanmean(timeDifferencesPerCycle,3);
        
        numCyclesCofiringThisDir=numCyclesCofiring(:,:,di);
        
        for u=1:numUnits
            noOverlapThisRow=numCyclesCofiringThisDir(u,:)<numCyclesThreshPerRow(u);
            numCyclesCofiringThisDir(u,:)=numCyclesCofiringThisDir(u,:)/ numCyclesMaxPerRow(u);
            numCyclesCofiringThisDir(u,noOverlapThisRow)=NaN;
            %numCyclesCofiringThisDir(noOverlapThisRow,u)=NaN;
            
            avgTimeDiffPerCycle(u,noOverlapThisRow)=NaN;
            %avgTimeDiffPerCycle(noOverlapThisRow,u)=NaN;
            
        end
        
        avgTimeDiffPerCycle=avgTimeDiffPerCycle(sortedMeanCenterIdxPerDi(:,di),sortedMeanCenterIdxPerDi(:,di));
        numCyclesCofiringThisDir=numCyclesCofiringThisDir(sortedMeanCenterIdxPerDi(:,di),sortedMeanCenterIdxPerDi(:,di));
        
        removeRowIdxes=[];
        removeColIdxes=[];
        for n=1:numUnits
            if(isnan(max(avgTimeDiffPerCycle(n,:))) || sum(~isnan(avgTimeDiffPerCycle(n,:)))==1)
                removeRowIdxes=[removeRowIdxes; n];
            end
            
            if(isnan(max(avgTimeDiffPerCycle(:,n))) || sum(~isnan(avgTimeDiffPerCycle(:,n)))==1)
                removeColIdxes=[removeColIdxes; n];
            end
        end
        
        
        avgTimeDiffPerCycle(removeRowIdxes,:)=[];
        avgTimeDiffPerCycle(:,removeColIdxes)=[];
        numCyclesCofiringThisDir(removeRowIdxes,:)=[];
        numCyclesCofiringThisDir(:,removeColIdxes)=[];
        
        figure; 
        subplot(1,2,1)
        imagesc(avgTimeDiffPerCycle)
        colormap(jet)
        colorbar
        caxis([-0.065 0.065])
        
        subplot(1,2,2)
        colormap(gca,parula)
        imagesc(numCyclesCofiringThisDir)
        colorbar
        %}
        
    %end   
    close all
end
        
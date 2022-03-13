close all; clear all; clc

processedDataDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData';
timeBoundSaveDir='./hc3Images/timeVsPhaseWithBound';

logStatsInfo=load('phaseVsTimeAndLogTimeStats.mat');
fieldIsUsedForLogTimeStats=logStatsInfo.fieldIsUsedForLogTimeStats;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%use all fields
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fieldIsUsedForLogTimeStats=ones(size(fieldIsUsedForLogTimeStats));

setTightSubplots


mycmap=getBiColorMap();

touchDir(timeBoundSaveDir)

unitDataDir=fullfile(processedDataDir,'unitSpikeInfoStructs');

unitFilePaths=getFilePathsRegex(unitDataDir,'*mat');


sortData=load('speedPhaseRsortData.mat');

speedPhaseRSortIdx=sortData.speedPhaseRSortIdx;

avgPerTraversalMinNumSpikes=1;
%avgPerTraversalMinNumSpikes=0;
minSpeedSTD=0.1;
%minSpeedSTD=0.2;

maxRunningSpeed=0.6;
maxRunningSpeed=0.75;
maxRunningSpeed=0.8;

plotSpeedInField=1;
plotSpeedInField=0;
if(plotSpeedInField)
    timeVsPhase=figure
    speedVsOffset=figure
end

singleCellLowSpeedMax=0.2;
singleCellHighSpeedMin=0.5;


%singleCellLowSpeedMax=0.15;
%singleCellHighSpeedMin=0.55;
singleCellLowSpeedMax=0.2;
singleCellHighSpeedMin=0.6;

singleCellLowSpeedMax=0.2;
singleCellHighSpeedMin=0.6;

singleCellLowSpeedMax=0.2;
singleCellHighSpeedMin=0.55;


singleCellLowSpeedMax=0.2;
singleCellHighSpeedMin=0.6;


%singleCellLowSpeedMax=0.2;
%singleCellHighSpeedMin=0.5;

singleCellLowSpeedMax=0.2;
singleCellHighSpeedMin=0.5;
singleCellHighSpeedMin=0.6;

singleCellLowSpeedMax=0.15;
%singleCellHighSpeedMin=0.5;
singleCellHighSpeedMin=0.7;

singleCellLowSpeedMax=0.3;
%singleCellHighSpeedMin=0.5;
singleCellHighSpeedMin=0.8;

singleCellLowSpeedMax=0.25;
%singleCellHighSpeedMin=0.5;
singleCellHighSpeedMin=0.7;

singleCellLowSpeedMax=0.25;
%singleCellHighSpeedMin=0.5;
singleCellHighSpeedMin=0.5;

singleCellLowSpeedMax=0.3;
%singleCellHighSpeedMin=0.5;
singleCellHighSpeedMin=0.5;


singleCellLowSpeedMax=0.2;
%singleCellHighSpeedMin=0.5;
singleCellMidSpeedMin=0.45;
singleCellHighSpeedMin=0.7;


singleCellLowSpeedMax=0.2;
%singleCellHighSpeedMin=0.5;
singleCellMidSpeedMin=0.4;
singleCellHighSpeedMin=0.6;

singleCellLowSpeedMax=0.15;
%singleCellHighSpeedMin=0.5;
singleCellMidSpeedMin=0.4;
singleCellHighSpeedMin=0.75;

singleCellLowSpeedMax=0.2;
%singleCellHighSpeedMin=0.5;
singleCellMidSpeedMin=0.4;
singleCellHighSpeedMin=0.6;
singleCellHighSpeedMax=0.7;

%speedCategoryBounds=(1:7)/10; %6 categories, 7 bounds
speedCategoryBounds=(2:7)/10; %5 categories, 6 bounds
speedCategoryBounds=(2.5:1:7.5)/10; %4 categories, 5 bounds

speedCategoryBounds=[0.2 0.4 0.6 0.8];

speedCategoryBounds=[0.25 0.35 0.55 0.65];

speedCategoryBounds=[0.15 0.25 0.35 0.45 0.55 0.65];

speedCategoryBounds=[0.35 0.45 0.55 0.65 0.75];

speedCategoryBounds=[0.15 0.25 0.35 0.45 0.55 0.65 0.75];

speedCategoryBounds=[0.2 0.3 0.4 0.5 0.6 0.7 0.8];

speedCategoryBounds=[0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7];

speedCategoryBounds=[0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7];

speedCategoryBounds=[0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7];

speedCategoryBounds=0.35:0.025:0.65; %high res (2.5cnormalized bin width )
speedCategoryBounds=0.35:0.02:0.65; %high res (2cnormalized bin width )
speedCategoryBounds=0.4:0.02:0.7; %high res (2cnormalized bin width )
speedCategoryBounds=0.35:0.01:0.7; %high res (2cnormalized bin width )
speedCategoryBounds=0.4:0.01:0.65; %high res (2cnormalized bin width )
speedCategoryBounds=0.3:0.02:0.65; %high res (2cnormalized bin width )
%speedCategoryBounds=[0 0.25 0.35 0.45 0.55 0.65 0.75];
speedCategoryBounds=0.35:0.01:0.7; %high res (2cnormalized bin width )
speedCategoryBounds=0.4:0.01:0.65; %high res (2cnormalized bin width )
speedCategoryBounds=0.4:0.01:0.6; %high res (2cnormalized bin width )
speedCategoryBounds=0.3:0.01:0.6; %high res (2cnormalized bin width )
speedCategoryBounds=0.35:0.025:0.7; %high res (2cnormalized bin width )
speedCategoryBounds=0.35:0.025:0.7; %high res (2cnormalized bin width )
%speedCategoryBounds=0.35:0.05:0.7; %high res (2cnormalized bin width )
%speedCategoryBounds=0.2:0.15:0.8;

%speedCategoryBounds=0.1:0.15:0.7; %4 categories, 5 bounds
   
speedCategoryBounds=0.01:0.05:0.25; %high res (1/100 to 1/4 fields/cycle bin width )
speedCategoryBounds=0.05:0.025:0.2; %high res (1/100 to 1/4 fields/cycle bin width )
speedCategoryBounds=0.08:0.02:0.16; %high res (1/100 to 1/4 fields/cycle bin width )

speedCategoryBounds=0.08:0.01:0.16; %high res (1/100 to 1/4 fields/cycle bin width )
speedCategoryBounds=0.04:0.02:0.18; %high res (1/100 to 1/4 fields/cycle bin width )
speedCategoryBounds=0.08:0.02:0.16; %high res (1/100 to 1/4 fields/cycle bin width )

speedCategoryBounds=1/100*[8 10 11 12 13 14 16];

speedCategoryBounds=1/100*[ 8 10 12 14 15 16 17 19.5 20.5 24 26];
speedCategoryBounds=1/100*[ 8 10 12 14 15 16 17 19.5 20.5];
speedCategoryBounds=1/100*[ 8 10 12 14 15 16 17 ];
speedCategoryBounds=1/100*[ 6 8 10 12 14 15 16 17 ];



speedCategoryBounds=0.35:0.025:0.7; %high res (2cnormalized bin width )
speedCategoryBounds=0.35:0.05:0.7; %high res (2cnormalized bin width )

speedCategoryBounds=0.45:0.05:0.75; %high res (2cnormalized bin width )
speedCategoryBounds=0.5:0.025:0.75; %high res (2cnormalized bin width )
speedCategoryBounds=0.45:0.05:0.7; %high res (2cnormalized bin width )
speedCategoryBounds=0.35:0.05:0.7; %high res (2cnormalized bin width )
%speedCategoryBounds=0.05:0.05:0.35; %high res (2cnormalized bin width )
speedCategoryBounds=1/100*[ 6 7 8 9 10 11 12  ];

%speedCategoryBounds=0.6:0.1:1.4;

speedCategoryBounds=0.7:0.1:1.3;

speedCategoryBounds=0.6:0.15:1.4;

speedCategoryBounds=0.6:0.15:1.4;
speedCategoryBounds=0.6:0.1:1.1;

speedCategoryBounds=0.6:0.1:1.3;
speedCategoryBounds=0.6:0.05:1;

speedCategoryBounds=0.65:0.1:1.25;
speedCategoryBounds=0.65:0.1:1.15;

%speedCategoryBounds=0.65:0.05:1;
speedCategoryBounds=0.65:0.05:0.95;

speedCategoryBounds=0.8:0.05:1.2;
speedCategoryBounds=0.8:0.05:1.2;

speedCategoryBounds=0.8:0.05:1.15;
speedCategoryBounds=0.8:0.075:1.2;

restrictToSameSizeFields=0;

%minFieldSizeThetaCycles=6;
%maxFieldSizeThetaCycles=14;

minFieldSizeThetaCycles=7;
maxFieldSizeThetaCycles=13;


minFieldSizeThetaCycles=6;
maxFieldSizeThetaCycles=14;

minFieldSizeThetaCycles=7;
maxFieldSizeThetaCycles=15;


minFieldSizeThetaCycles=6;
maxFieldSizeThetaCycles=20;



%minFieldSizeThetaCycles=6;
%maxFieldSizeThetaCycles=14;


%{
minFieldSizeThetaCycles=4;
maxFieldSizeThetaCycles=7;

minFieldSizeThetaCycles=15;
maxFieldSizeThetaCycles=20;
%}


%{
minFieldSizeThetaCycles=7;
maxFieldSizeThetaCycles=13;
%}

%{
speedCategoryBounds=0.3:0.1:2;
speedCategoryBounds=0.1:0.2:1.3;
speedCategoryBounds=0.6:0.1:1.2;
%}


numSpeedGroups=length(speedCategoryBounds)-1;


allPhaseVsTimePerSpeedGroupAllFields=cell(numSpeedGroups,1);


%{
singleCellLowSpeedMax=0.25;
singleCellHighSpeedMin=0.45;
%}


%{
singleCellLowSpeedMax=0.3;
singleCellHighSpeedMin=0.4;
%}

%{
singleCellLowSpeedMax=0.25;
singleCellHighSpeedMin=0.45;
%}

%{
singleCellLowSpeedMax=0.2;
singleCellHighSpeedMin=0.6;
%}


%maxRunningSpeed=0.5;
%maxRunningSpeed=0.45;
%maxRunningSpeed=0.8;
%maxRunningSpeed=0.5;
%maxRunningSpeed=0.7;
%maxRunningSpeed=0.65;
%maxRunningSpeed=0.55;
minRunningSpeed=0.05;
minRunningSpeed=0.1;
minRunningSpeed=0.2;
minRunningSpeed=0.25;
%minRunningSpeed=0.15;

minRunningSpeed=0.1;
%minRunningSpeed=0.075;

maxSpeedTimeCorr=0.1;
maxSpeedTimeCorr=0.05;
maxSpeedTimeCorr=1;
%maxSpeedTimeCorr=0.025;

speedPerSpike=1;
speedPerSpike=0;

plotIndFields=1;
plotIndFields=0;

totalFieldCount=0;
totalFieldCountWithSpeedVar=0;
%setTightSubplots
maximizeInitializePhase=0;
%maximizeInitializePhase=1;

numSpeedBins=60;
numSpeedBins=20;
%numSpeedBins=10;
numSpeedBins=15;
numSpeedBins=20;
numSpeedBins=10;
numSpeedBins=15;
%numTimeBins=40;
%numSpeedBins=5
%numSpeedBins=3
numSpeedBins=10;
numSpeedBins=15;
numSpeedBins=5;
numSpeedBins=7;
numSpeedBins=10;
numSpeedBins=5;

numSpeedBins=8;
numSpeedBins=12;
numSpeedBins=10;

minNumSpikes=500;

minNumSpikes=1000;

minNumPtsCorr=8;

%minNumSpikes=900;
minNumSpikes=0;
%minNumSpikes=500;
%minNumSpikes=1000;

minRunSpeedDisp=minRunningSpeed; %normalized
maxRunSpeedDisp=maxRunningSpeed; %normalized

plotByLap=0; 
plotBySpeed=0;
plotByDist=1;



%maxSpeedSEM=0.01;

maxSpeedRange=0.2;

maxSpeedRange=0.3;

maxSpeedRange=0.1;
maxSpeedRange=0.15;

maxSpeedRange=0.2;
maxSpeedRange=0.15;
maxSpeedRange=0.1;
%maxSpeedRange=0.25;
maxSpeedRange=0.15;
maxSpeedRange=0.125;

maxSpeedRange=0.1;

maxSpeedRange=0.05;
maxSpeedRange=0.025;
maxSpeedRange=0.03;
maxSpeedRange=0.025;
maxSpeedRange=0.03;
maxSpeedRange=0.05;
maxSpeedRange=0.1;
maxSpeedRange=0.025;


minSpeedStdvAcrossLaps=0.15; % 79 fields
minSpeedStdvAcrossLaps=0.08; % 798 fields
minSpeedStdvAcrossLaps=0.1; % ~489 fields
minSpeedStdvAcrossLaps=0; % 1227 fields

allUnitSpikeSpeeds=[];
allUnitSpikePhases=[];
allLowSpeedPhaseVsTimeFracs=[];
allMidSpeedPhaseVsTimeFracs=[];
allHighSpeedPhaseVsTimeFracs=[];
allHighestSpeedPhaseVsTimeFracs=[];

allUnitSpikeTimesInField=[];
allSpeedPhaseR=[];
allSpeedPhaseSlopes=[];
allRangesPerTraversal=[];

allHighVsLow=[];
allCellPairedLowHigh=[];

figure

 spikePhasesPerField={};
 spikeTimeFracInFieldPerField={};
 spikeTimesInExpPerField={};
 unitInfoPathPerField={};
 dirPerField={};

startFi=1;
allIndivFieldRs=[];
allIndivFieldSlopes=[];

includedFieldCount=0;


for fi=startFi:length(unitFilePaths)
   
    if(mod(fi,50)==0)
        disp(fi)
    end
    currUnitFilePath=unitFilePaths{fi};
    currUnitStruct=load(currUnitFilePath);
    
    if(isfield(currUnitStruct,'cycleTimeFieldBoundPerDirPerField'))
        currUnitData=currUnitStruct.unitInfo;

        currUnitFileName=getFileNameFromPath(currUnitFilePath);
        
        cycleTimeFieldBoundPerDirPerField=currUnitStruct.cycleTimeFieldBoundPerDirPerField;

        numRightwardFields=length(currUnitStruct.rightwardFieldStartEndM)/2;
        numLeftwardFields=length(currUnitStruct.leftwardFieldStartEndM)/2;

        maxNumFieldsPerDir=max([numRightwardFields numLeftwardFields]);

        %if(redoBounds)
            rightwardTimeBounds=NaN(maxNumFieldsPerDir,1);
        %end
        
        for di=1:2
            if(di==1)
                currFieldDirStr='rightward';
                numFields=numRightwardFields;
                
            else
                currFieldDirStr='leftward';
                numFields=numLeftwardFields;
            end
            
            for fii=1:numFields
                currFieldThetaData=currUnitStruct.thetaBasedTimeVsPhaseInfo.(currFieldDirStr).(sprintf('field%d',fii));

                allSpikesGoodLapIdxes=currFieldThetaData.allSpikeGoodLapIdxes;
                
                if(length(allSpikesGoodLapIdxes)<minNumSpikes)
                    continue
                end
                
                %currFieldTimeBound=currUnitStruct.cycleTimeFieldBoundPerDirPerField{fii,di}-1; %-1 becaues manually entered before changing first theta cycle from 1 to 0
                
                %currFieldTimeBound=currFieldTimeBound*0.9*0.9*0.95;
                
                if(speedPerSpike)
                    allSpikesCurrSpeeds=currFieldThetaData.allLapInFieldSpikeSpeeds;
                else
                    allSpikesCurrSpeeds=abs(currFieldThetaData.allLapInFieldSpikeTraversalAvgSpeeds);
                end

                allSpikesCurrTimeFieldPhases=currFieldThetaData.allLapInFieldSpikePhases;
                
                allSpikeExpTimes=currFieldThetaData.allLapInFieldSpikeTimesInExp;
                
                meanSpeedPerTraversal=currFieldThetaData.allLapInFieldMeanSpeedPerTraversal;
                meanPhasePerTraversal=currFieldThetaData.allLapInFieldMeanPhasePerTraversal;
                numCyclesPerTraversal=currFieldThetaData.numThetaCyclesInFieldPerLap;
                %currFieldTimeBound=(currUnitStruct.cycleTimeFieldBoundPerDirPerField{fii,di}-1)*0.8; %-1 becaues manually entered before changing first theta cycle from 1 to 0
                currFieldTimeBound=(currUnitStruct.cycleTimeFieldBoundPerDirPerField{fii,di}-1); %-1 becaues manually entered before changing first theta cycle from 1 to 0
                        %currFieldTimeBound=(currUnitStruct.cycleTimeFieldBoundPerDirPerField{fii,di}); %-1 becaues manually entered before changing first theta cycle from 1 to 0

                
                        if(restrictToSameSizeFields)
                            if(currFieldTimeBound>maxFieldSizeThetaCycles || currFieldTimeBound<minFieldSizeThetaCycles)
                                continue
                            end
                        end
                        
                if(std(meanSpeedPerTraversal)<minSpeedStdvAcrossLaps)
                    continue
                end
                
                semSpeedPerTraversal=currFieldThetaData.allLapInFieldSpeedRangePerTraversal;
                
                allRangesPerTraversal=[allRangesPerTraversal;semSpeedPerTraversal(:)];
                %meanRangePerTraversal=
                
                 uniformSpeedPerTraversal=semSpeedPerTraversal<=maxSpeedRange;
                 
   
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 %TAKE ONLY LAPS WITH ~UNIFORM SPEED THROUGH FIELD
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 %meanSpeedPerTraversal=meanSpeedPerTraversal(uniformSpeedPerTraversal);
                 %meanSpeedPerTraversal=1./numCyclesPerTraversal(uniformSpeedPerTraversal);
                 
                 meanSpeedPerTraversal=currFieldTimeBound./numCyclesPerTraversal(uniformSpeedPerTraversal);
                 
                  meanPhasePerTraversal=meanPhasePerTraversal(uniformSpeedPerTraversal);
                 %{
                  plot(meanSpeedPerTraversal,meanPhasePerTraversal,'k.');
                 drawnow
                 hold on
                  %}
                 
                
                 
                 
                
                 
                 %{
                  lowSpeedLapIdxes=meanSpeedPerTraversal<=singleCellLowSpeedMax;
                  midSpeedLapIdxes= meanSpeedPerTraversal>singleCellLowSpeedMax &  meanSpeedPerTraversal<=singleCellMidSpeedMin;
                  highSpeedLapIdxes=meanSpeedPerTraversal>singleCellMidSpeedMin & meanSpeedPerTraversal<=singleCellHighSpeedMin;
                  highestSpeedLapIdxes=meanSpeedPerTraversal>singleCellHighSpeedMin;
                 
                  lowSpeedOnlyAllPhaseVsTime=vertcat(currFieldThetaData.allLapInFieldAllTimeVsPhasePerTraversal{lowSpeedLapIdxes});
                  midSpeedOnlyAllPhaseVsTime=vertcat(currFieldThetaData.allLapInFieldAllTimeVsPhasePerTraversal{midSpeedLapIdxes});
                  highSpeedOnlyAllPhaseVsTime=vertcat(currFieldThetaData.allLapInFieldAllTimeVsPhasePerTraversal{highSpeedLapIdxes});
                  highestSpeedOnlyAllPhaseVsTime=vertcat(currFieldThetaData.allLapInFieldAllTimeVsPhasePerTraversal{highestSpeedLapIdxes});

                 %}
                  
                 allPhaseVsTimePerSpeedGroupCurrField=cell(numSpeedGroups,1);
                 currSpeedGpIdxes=cell(numSpeedGroups,1);
                    
                  for si=1:numSpeedGroups
                      currSpeedBinMin=speedCategoryBounds(si);
                      currSpeedBinMax=speedCategoryBounds(si+1);
                      
                      currSpeedGpIdxes{si}=meanSpeedPerTraversal>=currSpeedBinMin & meanSpeedPerTraversal<currSpeedBinMax;
                      
                      allPhaseVsTimePerSpeedGroupCurrField{si}=vertcat(currFieldThetaData.allLapInFieldAllTimeVsPhasePerTraversal{currSpeedGpIdxes{si}});
                  end
            
              

                currFieldTimeBound=currUnitStruct.cycleTimeFieldBoundPerDirPerField{fii,di}-1; %-1 becaues manually entered before changing first theta cycle from 1 to 0
                
                currFieldTimeBound=currFieldTimeBound*0.9*0.9*0.95;
                
                allSpikesCurrTimeFieldFracs=currFieldThetaData.allLapSpikeCycleIDsFromFieldEntry/currFieldTimeBound;
                        
                goodLapSpikesCurrTimeFieldFracs=allSpikesCurrTimeFieldFracs(allSpikesGoodLapIdxes);
                goodLapSpikesCurrSpeeds=allSpikesCurrSpeeds(allSpikesGoodLapIdxes);
                goodLapSpikesCurrTimeFieldPhases=allSpikesCurrTimeFieldPhases(allSpikesGoodLapIdxes);
                
                allSpikeExpTimes=allSpikeExpTimes(allSpikesGoodLapIdxes);
              
                allPhases=goodLapSpikesCurrTimeFieldPhases(:);
               

                %{
                plot(goodLapSpikesCurrTimeFieldFracs,allPhases,'k.')
                hold on
                xlim([0 1])
                ylim([0 360])
                %}
                
                if(plotIndFields & fieldIsUsedForLogTimeStats(totalFieldCount+1))
                    subplot(15,12,includedFieldCount+1)
                    %{
                    plot(goodLapSpikesCurrSpeeds,allPhases,'k.','MarkerSize',1)
                    xlabel('Time in field (in cycles, frac.)')
                    
                    
                    ylabel('Theta phase (deg)')
                    %}
                    
                    plot(meanSpeedPerTraversal,meanPhasePerTraversal,'k.','MarkerSize',10)
                    xlabel('Mean speed in traversal (normalized)')
                    
                    
                    ylabel('Mean theta phase in traversal (deg)')
                    xlim([0 1])
                    ylim([0 360])
                    drawnow
                    
                   
                end
                
                if(fieldIsUsedForLogTimeStats(totalFieldCount+1))
                     notNaNIdxes=~isnan(meanSpeedPerTraversal) & ~isnan(meanPhasePerTraversal);
                    [m,b,R]=getLinearFit(meanSpeedPerTraversal(notNaNIdxes), meanPhasePerTraversal(notNaNIdxes));
                    includedFieldCount=includedFieldCount+1;
                    allIndivFieldRs=[allIndivFieldRs; R];
                    allIndivFieldSlopes=[allIndivFieldSlopes; m];
                end
                
                corrIdxes=(goodLapSpikesCurrTimeFieldFracs<1);
                %{
                rMat=corrcoef(goodLapSpikesCurrSpeeds(corrIdxes),goodLapSpikesCurrTimeFieldFracs(corrIdxes));
                
                if(length(rMat)==1) %not enough points
                    continue
                end
                speedTimeCorr=rMat(2,1);
                
                if(abs(speedTimeCorr)>maxSpeedTimeCorr)
                    continue
                end
                %}
                if(min(goodLapSpikesCurrTimeFieldFracs)>maxSpeedTimeCorr)
                    continue
                end
                
                if(avgPerTraversalMinNumSpikes)
                    if(fieldIsUsedForLogTimeStats(totalFieldCount+1))
                        allUnitSpikeSpeeds=[allUnitSpikeSpeeds; meanSpeedPerTraversal];
                        allUnitSpikePhases=[allUnitSpikePhases; meanPhasePerTraversal];
                        
                        %{
                        allLowSpeedPhaseVsTimeFracs=[allLowSpeedPhaseVsTimeFracs;lowSpeedOnlyAllPhaseVsTime];
                        allMidSpeedPhaseVsTimeFracs=[allMidSpeedPhaseVsTimeFracs;midSpeedOnlyAllPhaseVsTime];
                        allHighSpeedPhaseVsTimeFracs=[allHighSpeedPhaseVsTimeFracs;highSpeedOnlyAllPhaseVsTime];
                        allHighestSpeedPhaseVsTimeFracs=[allHighestSpeedPhaseVsTimeFracs;highestSpeedOnlyAllPhaseVsTime];
                        %}
                        for si=1:numSpeedGroups
                            allPhaseVsTimePerSpeedGroupAllFields{si} = [ allPhaseVsTimePerSpeedGroupAllFields{si}; allPhaseVsTimePerSpeedGroupCurrField{si}];
                        end
                        
                        % allUnitSpikeTimeInField=[allUnitSpikeTimeInField; meanTimePerTraversal];
                    end 
                     %allUnitSpikeTimesInField=[allUnitSpikeTimesInField; ];
                else
                    if(fieldIsUsedForLogTimeStats(totalFieldCount+1))
                        allUnitSpikeSpeeds=[allUnitSpikeSpeeds; goodLapSpikesCurrSpeeds(:)];
                        allUnitSpikePhases=[allUnitSpikePhases; allPhases(:)];
                         allUnitSpikeTimesInField=[allUnitSpikeTimesInField; goodLapSpikesCurrTimeFieldFracs(:)];

                         spikePhasesPerField{totalFieldCount+1}=allPhases(:);
                         spikeTimeFracInFieldPerField{totalFieldCount+1}=goodLapSpikesCurrTimeFieldFracs(:);
                         spikeTimesInExpPerField{totalFieldCount+1}=allSpikeExpTimes(:);


                         if(length(allPhases) ~= length(allSpikeExpTimes))
                             disp('')

                         end
                         unitInfoPathPerField{totalFieldCount+1}=currUnitFilePath;
                         dirPerField{totalFieldCount+1}=di;
                     
                    end
                end
                
                totalFieldCount=totalFieldCount+1;
                
                hasSpeedVariation=nanstd(meanSpeedPerTraversal)>=minSpeedSTD & (sum(~isnan(meanSpeedPerTraversal) & ~isnan(meanPhasePerTraversal))>minNumPtsCorr);
                
                 %speedPhaseR=getCorrCoeff(meanSpeedPerTraversal,meanPhasePerTraversal);
                 [speedPhaseR,p,slopeDegPerSpeed]=getCircCorrCoeff(meanSpeedPerTraversal,meanPhasePerTraversal);
                 
                   minNumLaps=1;
                   minNumLaps=3;
                   %minNumLaps=5;
                   %minNumLaps=4;
                   %minNumLaps=2;
                   minNumLaps=5;
                   minNumLaps=1;
                       %minNumLaps=2;
                  
                      
                       
                       lowSpeedLapIdxes=currSpeedGpIdxes{1};
                       highSpeedLapIdxes=currSpeedGpIdxes{end};
                       
                  if(sum(lowSpeedLapIdxes)>=minNumLaps && sum(highSpeedLapIdxes)>=minNumLaps)
                   %if(sum(currSpeedGpIdxes{1})>=minNumLaps && sum(currSpeedGpIdxes{end})>=minNumLaps)

                  
                      %lowSpeedLapAvgPhase=circMeanDeg(meanPhasePerTraversal(lowSpeedLapIdxes));
                     % highSpeedLapAvgPhase=circMeanDeg(meanPhasePerTraversal(highSpeedLapIdxes));
                      
                      lowSpeedLapAvgPhase=nanmean(meanPhasePerTraversal(lowSpeedLapIdxes));
                      highSpeedLapAvgPhase=nanmean(meanPhasePerTraversal(highSpeedLapIdxes));

                      %highVsLowCircDiff=angdiffDeg([lowSpeedLapAvgPhase highSpeedLapAvgPhase ]);

                       highVsLowDiff=diff([lowSpeedLapAvgPhase highSpeedLapAvgPhase ]);

                      cellPairedLowHigh=[lowSpeedLapAvgPhase highSpeedLapAvgPhase ];
                      
                      if(length(cellPairedLowHigh)==2)
                      
                          allCellPairedLowHigh=[allCellPairedLowHigh; cellPairedLowHigh];
                          %allHighVsLow=[allHighVsLow;highVsLowCircDiff];
                          allHighVsLow=[allHighVsLow;highVsLowDiff];
                      end
                  end
                 
                 if(hasSpeedVariation)
                     allSpeedPhaseR=[allSpeedPhaseR; speedPhaseR];
                     allSpeedPhaseSlopes=[allSpeedPhaseSlopes; slopeDegPerSpeed];
                     totalFieldCountWithSpeedVar=totalFieldCountWithSpeedVar+1;
                    goodVarSpeedPhaseR(totalFieldCountWithSpeedVar)=speedPhaseR;
                 end
                
                 
                 cmapLength=64;
                 slopeColors=jet(cmapLength);
                 %cmapLength=64;
                 %slopeColors=mycmap; 
                if(plotSpeedInField && hasSpeedVariation)
                    
                   
                    figure(speedVsOffset)
                    plotIdx=find(speedPhaseRSortIdx==totalFieldCountWithSpeedVar);
                    %subplot(20,25,totalFieldCountWithSpeedVar)
                    %subplot(20,25,plotIdx)
                    subplot(18,25,plotIdx)
                    
                    %plot(goodLapSpikesCurrTimeFieldFracs,goodLapSpikesCurrSpeeds,'k.','MarkerSize',2)
                   
                    %colorIdx=min(cmapLength,round((slopeDegPerSpeed+1000)/2000*cmapLength));
                    colorIdx=min(cmapLength,round((1+speedPhaseR)/2*cmapLength));

                    
                    if(colorIdx<1)
                        colorIdx=1;
                    end
                    currSlopeColor=slopeColors(colorIdx,:);
                    
                    figure(speedVsOffset)
                    %plot(meanSpeedPerTraversal,meanPhasePerTraversal,'k.','MarkerSize',2)
                    plot(meanSpeedPerTraversal,meanPhasePerTraversal,'.','Color',currSlopeColor,'MarkerSize',8)
                   
                    box off
                    %xlabel('Time in field (in cycles, frac.)')
                    %ylabel('Speed (normalized)')
                    %xlim([0 1])
                    %ylim([0 0.8])
                    xlim([0 0.8])
                    ylim([0 360])
                    
                    if(mod(plotIdx,25)~=1)
                         yticklabels({})
                    else
                        ylabel('\mu Phase (\circ)')
                    end
                    if(plotIdx<(18*25-24))
                        xticklabels({})
                    else
                        xlabel('Speed (normalized)')
                    end
                    
                    %set(gca,'Color','k')
                    maxFig
                     setFigFontTo(9)
                    drawnow
                    
                    figure(timeVsPhase)
                    subplot(18,25,plotIdx)
                     plot(goodLapSpikesCurrTimeFieldFracs,allPhases,'.','Color',currSlopeColor,'MarkerSize',5)
                    box off
                    %xlabel('Time in field (in cycles, frac.)')
                    %ylabel('Speed (normalized)')
                    %xlim([0 1])
                    %ylim([0 0.8])
                    xlim([0 1])
                    ylim([0 360])
                    
                    if(mod(plotIdx,25)~=1)
                         yticklabels({})
                    else
                        ylabel('\Theta Phase (\circ)')
                    end
                    if(plotIdx<(18*25-24))
                        xticklabels({})
                    else
                        xlabel('Time (frac)')
                    end
                    
                    
                    maxFig
                    setFigFontTo(9)
                    drawnow
                    
                end
               
                
                %totalFieldCount=totalFieldCount+1;
                 
            
            end
        end
        
    end
    
end

%save('allUnitInFieldGoodSpeedSpikePhasesAndTimes.mat','allUnitSpikeTimesInField','allUnitSpikePhases','totalFieldCount','spikePhasesPerField','spikeTimeFracInFieldPerField','spikeTimesInExpPerField','unitInfoPathPerField','dirPerField')
disp('saved.')

%addAutomatedTimeFieldBoundNormTimes

includedFieldCount


%%
%{
onlyLateField=1;
if(onlyLateField)
    %lateSpikeIdxes=allUnitSpikeTimesInField<0.5;
     %lateSpikeIdxes=allUnitSpikeTimesInField<0.6;
      lateSpikeIdxes=allUnitSpikeTimesInField>0.6;
      lateSpikeIdxes=allUnitSpikeTimesInField>0;
    allUnitSpikeSpeeds(~lateSpikeIdxes)=NaN;
    allUnitSpikePhases(~lateSpikeIdxes)=NaN;
end
%}



%{
fG=figure

numTimeBins=30;
timeEdges=linspace(0,1,numTimeBins+1);
timeBinCenters=edgesToBins(timeEdges);

[pxySmoothSpeedVsTime]=getJointDistr(allUnitSpikeSpeeds,allUnitSpikeTimesInField,speedEdges,timeEdges,fG);


figure;
timeBinColors=jet(numSpeedBins)
for ti=1:numSpeedBins
    currPdist=circSmooth(pxySmoothSpeedVsTime(ti,:),5);
    currPdist=currPdist/sum(currPdist);
    plot(timeBinCenters,currPdist,'Color',timeBinColors(ti,:),'LineWidth',3)
    
    %[~,maxID]=max(currPdist);
    %maxProbPhasePerTimeBin(ti)=phaseBinCenters(maxID);
    

    %circMeanPhasePerSpeedBin(ti)=circMeanDeg(
    
    
    hold on
end
colormap(gca,jet)
cb=colorbar
ylabel(cb,'Speed (normalized)')
xlabel('Time in field (frac)')
ylabel('Probability')
caxis([minRunningSpeed maxRunningSpeed])
xlim([0 1])

title('Running speed vs time in field distribution')
maxFig
setFigFontTo(18)
saveas(gcf,'speedInFieldVsThetaPhaseDistr.tif')
%print('-r50',sprintf('speedInFieldVsTimeInFieldDistr'),'-dpng')




%}
%%
fH=figure;

%{
numTimeBins=20;
numTimeBins=40;
numTimeBins=60;
numTimeBins=120;
numTimeBins=20;
%}
%numSpeedBins=20
%numSpeedBins=7
%numSpeedBins=5
%numSpeedBins=10

minRunningSpeed=min(speedCategoryBounds);
maxRunningSpeed=max(speedCategoryBounds);
speedEdges=linspace(minRunningSpeed,maxRunningSpeed,numSpeedBins+1);
speedBinCenters=edgesToBins(speedEdges);
speedBinWidth=median(diff(speedBinCenters));

numPhaseBins=40;
numPhaseBins=60;
numPhaseBins=120;
numPhaseBins=80;
%numPhaseBins=20;
numPhaseBins=24;
numPhaseBins=30;
numPhaseBins=24;
numPhaseBins=60;
numPhaseBins=30;
numPhaseBins=20;
phaseEdges=linspace(0,360,numPhaseBins+1);
phaseBinCenters=edgesToBins(phaseEdges);

%[pxySmooth,pxy]=getJointDistr(allUnitSpikeSpeeds,allUnitSpikePhases,speedEdges,phaseEdges,fH);
[pxySmooth,pxy]=getJointDistrGivenX(allUnitSpikeSpeeds,allUnitSpikePhases,speedEdges,phaseEdges,fH);


title(sprintf('n=%d fields',totalFieldCount))
xlabel('Avg running speed in field traversal (normalized)')

if(avgPerTraversalMinNumSpikes)
    ylabel('Circ mean theta phase in field traversal (deg)')
else
    ylabel('Theta phase (deg)')
end


daspect([1 360 1])
%caxis([0 prctile(pxySmooth(:),97.5)])
%maxFig
setFigFontTo(18)
saveas(gcf,'speedVsThetaPhaseJointDistr.tif')


figure;
timeBinColors=jet(numSpeedBins)
maxProbPhasePerTimeBin=NaN(numSpeedBins,1);
for ti=1:numSpeedBins
    
    %{
    if(~(ti==1 || ti==numSpeedBins))
        continue
    end
    %}
    %currPdist=circSmooth(pxySmooth(ti,:),5);
    currPdist=circSmooth(pxy(ti,:),5);
    
    currPdist=currPdist/sum(currPdist);
    plot(phaseBinCenters,currPdist,'Color',timeBinColors(ti,:),'LineWidth',3)
    
    [~,maxID]=max(currPdist);
    maxProbPhasePerTimeBin(ti)=phaseBinCenters(maxID);
    

    %circMeanPhasePerSpeedBin(ti)=circMeanDeg(
    
    
    hold on
    drawnow
end
colormap(gca,jet)
cb=colorbar
ylabel(cb,'Avg running speed in field traversal (normalized)')

if(avgPerTraversalMinNumSpikes)
    xlabel('Circ mean theta phase in field traversal (deg)')
else
    xlabel('Theta phase (deg)')
end

ylabel('Probability')
xlim([0 360])
caxis([minRunningSpeed maxRunningSpeed])

title(sprintf('Avg running speed in traversal vs theta phase distribution, n=%d fields',totalFieldCount))
maxFig
setFigFontTo(18)
saveas(gcf,'speedInFieldVsThetaPhaseDistr.tif')
%print('-r50',sprintf('speedInFieldVsThetaPhaseDistr'),'-dpng')

figure;

plot(speedBinCenters,smooth(maxProbPhasePerTimeBin),'ko')
xlabel('Avg running speed in field traversal (normalized)')
ylabel('Most probable avg theta phase (deg)')
xlim([minRunningSpeed maxRunningSpeed])
ylim([0 360])

title(sprintf('n=%d fields',totalFieldCount))

daspect([1 360 1])
%maxFig
setFigFontTo(18)
maxFig
saveas(gcf,'SpeedVsMaxProbThetaPhase.tif')
%print('-r100',sprintf('SpeedVsMaxProbThetaPhase'),'-dpng')

%autoArrangeFigures()


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get circ mean phase per time bin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%numSpeedBinsForAvg=15;
%numSpeedBinsForAvg=7;

numSpeedBinsForAvg=numSpeedBins;



speedAvgEdges=linspace(minRunningSpeed,maxRunningSpeed,numSpeedBinsForAvg+1);
speedAvgBinCenters=edgesToBins(speedAvgEdges);
speedAvgBinWidth=median(diff(speedAvgBinCenters));

phaseMeanPerSpeedBin=NaN(numSpeedBinsForAvg,1);

for ti=1:numSpeedBinsForAvg
    currSpeedBinStart=speedAvgEdges(ti)-speedAvgBinWidth/2;
    currSpeedBinEnd=speedAvgEdges(ti)+speedAvgBinWidth/2;
    %currTimeBinCenter=timeBinCenters(ti);
    
    currBinIdxes=allUnitSpikeSpeeds>currSpeedBinStart & allUnitSpikeSpeeds<=currSpeedBinEnd;
    
    currBinPhases=allUnitSpikePhases(currBinIdxes);
    currBinMeanPhase=circMeanDeg(currBinPhases);
        %currBinMeanPhase=nanmean(currBinPhases);

    %currBinMeanPhase=circMedianDeg(currBinPhases);
    
    phaseMeanPerSpeedBin(ti)=currBinMeanPhase;
    
end

figure;
plot(speedAvgBinCenters,phaseMeanPerSpeedBin(1)-phaseMeanPerSpeedBin,'k.','MarkerSize',40)

%plot(speedAvgBinCenters,maxProbPhasePerTimeBin(1)-maxProbPhasePerTimeBin,'k.','MarkerSize',40)


%plot(speedAvgBinCenters,phaseMeanPerSpeedBin,'k.','MarkerSize',40)


xlabel('Avg running speed in field traversal (normalized)')
ylabel('Avg circ mean theta phase in field traversal (deg)')
xlim([minRunningSpeed maxRunningSpeed])
ylim([0 360])
ylim([0 180])
%ylim([100 300])
hold on; 
pL=plot(speedAvgBinCenters, getLogTime(speedAvgBinCenters,1.33333)*135,'r-','LineWidth',3)
%hold on; plot(speedAvgBinCenters, getLogTime(speedAvgBinCenters,1.33333)*125)
%hold on; plot(speedAvgBinCenters, getLogTime(speedAvgBinCenters,1.33333)*90)
%title(sprintf('Avg theta phase in field vs running speed, n=%d fields',totalFieldCount))
legend(pL,'log_{1.33}(speed)')
legend boxoff
title(sprintf('Average theta phase offset vs running speed, n=%d fields',sum(fieldIsUsedForLogTimeStats)))
setFigFontTo(18)
%maxFig

saveas(gcf,'speedVsCircMeanThetaPhase.tif')
%print('-r75',sprintf('speedVsCircMeanThetaPhase'),'-dpng')

%figure; histogram(allHighVsLow,20); 

%%
setTightSubplots_SpaceTime
timeEdges=linspace(0,1,31);
timeEdges=linspace(0,1,15);
timeEdges=linspace(0,1,31);
phaseEdges=linspace(0,360,31);
allPhasesHiVsLowSpeedFig=figure;

numUsedSpeedGroupsTotal=length(speedCategoryBounds)-1;
numUsedSpeedGroupsSoFar=0;
for si=1:numSpeedGroups
  currSpeedBinMin=speedCategoryBounds(si);
  currSpeedBinMax=speedCategoryBounds(si+1);
                      

  currSpeedGroupPhaseVsTime=allPhaseVsTimePerSpeedGroupAllFields{si};
  
  if(isempty(currSpeedGroupPhaseVsTime))
      continue
  end
  numUsedSpeedGroupsSoFar=numUsedSpeedGroupsSoFar+1;
    subplot(1,numUsedSpeedGroupsTotal,numUsedSpeedGroupsSoFar)
    withSmooth=1;
    dispNoise=normrnd(0,0.02,size(currSpeedGroupPhaseVsTime(:,1)));
    %try
    getJointDistr(currSpeedGroupPhaseVsTime(:,1)+dispNoise,currSpeedGroupPhaseVsTime(:,2),timeEdges,phaseEdges,allPhasesHiVsLowSpeedFig,withSmooth);
   % end
    %title(sprintf('speed < %.2f normalized',singleCellLowSpeedMax))
    %title(sprintf('%.2f field/cycle < speed < %.2f field/cycle',currSpeedBinMin,currSpeedBinMax))
        title(sprintf('%.2f normalized < speed < %.2f normalized',currSpeedBinMin,currSpeedBinMax))

    if(numUsedSpeedGroupsSoFar==1)
        climits=caxis;
    else
        caxis([climits])
    end

end

%{
subplot(1,4,2)
getJointDistr(allMidSpeedPhaseVsTimeFracs(:,1),allMidSpeedPhaseVsTimeFracs(:,2),timeEdges,phaseEdges,allPhasesHiVsLowSpeedFig);
title(sprintf('%.2f normalized < speed < %.2f normalized',singleCellLowSpeedMax,singleCellMidSpeedMin))

subplot(1,4,3)
getJointDistr(allHighSpeedPhaseVsTimeFracs(:,1),allHighSpeedPhaseVsTimeFracs(:,2),timeEdges,phaseEdges,allPhasesHiVsLowSpeedFig);
title(sprintf('%.2f normalized < speed < %.2f normalized',singleCellMidSpeedMin,singleCellHighSpeedMin))

subplot(1,4,4)
getJointDistr(allHighestSpeedPhaseVsTimeFracs(:,1),allHighestSpeedPhaseVsTimeFracs(:,2),timeEdges,phaseEdges,allPhasesHiVsLowSpeedFig);
title(sprintf('speed > %.2f normalized',singleCellHighSpeedMin))
%}
setFigFontTo(12)
%%
figure
%timeEdges=linspace(0,1,12);
timeEdges=linspace(0,1,8);
timeEdges=linspace(0,0.9,10);
timeEdges=linspace(0,1,7);
%timeEdges=linspace(0,0.9,8);
%timeEdges=linspace(0,1,6);
%timeEdges=linspace(0,1,12);

%speedGroupColors=copper(numSpeedGroups);

speedGroupColors=getBiColorMap(numSpeedGroups);

minCircShiftPhase=20;
for si=1:numSpeedGroups
%for si=2:numSpeedGroups
  currSpeedBinMin=speedCategoryBounds(si);
  currSpeedBinMax=speedCategoryBounds(si+1);
                      

  meanSpeed=mean([currSpeedBinMax currSpeedBinMin]);
  currSpeedGroupPhaseVsTime=allPhaseVsTimePerSpeedGroupAllFields{si};
  if(isempty(currSpeedGroupPhaseVsTime))
      continue
  end
  
  %speedColorIdx=round(numSpeedGroups*meanSpeed/max(speedCategoryBounds))
   speedColorIdx=si;


  
  [timeBinCenters,currSpeedGpPhaseAvgValues] = getBinnedCircAverages(currSpeedGroupPhaseVsTime(:,1),currSpeedGroupPhaseVsTime(:,2),timeEdges);
  
    currPhases=currSpeedGpPhaseAvgValues;
  
  shiftIdxes=currPhases<minCircShiftPhase;
  
  currPhases(shiftIdxes)=currPhases(shiftIdxes)+360;
    plot(timeBinCenters,currPhases,'-','Color', speedGroupColors(speedColorIdx,:),'LineWidth',5)
    hold on
        plot(timeBinCenters,currPhases,'.','Color', speedGroupColors(speedColorIdx,:),'MarkerSize',40)


end
%colormap(copper)
colormap(getBiColorMap)
cb=colorbar('north')

%ylabel(cb,'running speed (field fraction / cycle)')
ylabel(cb,'running speed (normalized)')
caxis([min(speedCategoryBounds) max(speedCategoryBounds)])
 xlabel('Time in field (frac)')
 ylabel('Mean theta phase (deg)')
 box off
 setFigFontTo(16)
 
 ylim([90 340])

%{
[timeBinCenters,lowSpeedPhaseAvgValues] = getBinnedCircAverages(allLowSpeedPhaseVsTimeFracs(:,1),allLowSpeedPhaseVsTimeFracs(:,2),timeEdges);
[timeBinCenters,midSpeedPhaseAvgValues] = getBinnedCircAverages(allMidSpeedPhaseVsTimeFracs(:,1),allMidSpeedPhaseVsTimeFracs(:,2),timeEdges);
[timeBinCenters,highSpeedPhaseAvgValues] = getBinnedCircAverages(allHighSpeedPhaseVsTimeFracs(:,1),allHighSpeedPhaseVsTimeFracs(:,2),timeEdges);
[timeBinCenters,highestSpeedPhaseAvgValues] = getBinnedCircAverages(allHighestSpeedPhaseVsTimeFracs(:,1),allHighestSpeedPhaseVsTimeFracs(:,2),timeEdges);

%subplot(1,4,1)
pL=plot(timeBinCenters,lowSpeedPhaseAvgValues,'k-','LineWidth',3)
hold on
plot(timeBinCenters,lowSpeedPhaseAvgValues,'ko');


%subplot(1,4,2)
hold on
pM=plot(timeBinCenters,midSpeedPhaseAvgValues,'b-','LineWidth',3)
hold on
plot(timeBinCenters,midSpeedPhaseAvgValues,'bo');

%subplot(1,4,3)
pH=plot(timeBinCenters,highSpeedPhaseAvgValues,'g-','LineWidth',3)

plot(timeBinCenters,highSpeedPhaseAvgValues,'go');


%subplot(1,4,4)
pHST=plot(timeBinCenters,highestSpeedPhaseAvgValues,'r-','LineWidth',3)

plot(timeBinCenters,highestSpeedPhaseAvgValues,'ro');
 ylim([50 360])
%}


%legend([pL,pM,pH,pHST],{sprintf('speed < %.2f normalized',singleCellLowSpeedMax), sprintf('%.2f normalized < speed < %.2f normalized',singleCellLowSpeedMax,singleCellMidSpeedMin),sprintf('%.2f normalized < speed < %.2f normalized',singleCellMidSpeedMin,singleCellHighSpeedMin),sprintf('speed > %.2f normalized',singleCellHighSpeedMin)})
%legend([pL,pM,pH,pHST],{sprintf('speed < %.2f normalized',singleCellLowSpeedMax), sprintf('%.2f normalized < speed < %.2f normalized',singleCellLowSpeedMax,singleCellMidSpeedMin),sprintf('%.2f normalized < speed < %.2f normalized',singleCellMidSpeedMin,singleCellHighSpeedMin),sprintf('speed > %.2f normalized',singleCellHighSpeedMin)})

%legend boxoff

%legend([p1 p2 p3 p4], {sprintf('speed < %.2f normalized',singleCellLowSpeedMax),...
%    sprintf('%.2f normalized < speed < %.2f normalized',singleCellLowSpeedMax,singleCellMidSpeedMin),sprintf('%.2f normalized < speed < %.2f normalized',singleCellMidSpeedMin,singleCellHighSpeedMin), sprintf('speed > %.2f normalized',singleCellHighSpeedMin)})



%save('speedVsOffsetData_Oct8_2021.mat')
%save('speedVsOffsetData_CloseTimeBound')
%save('speedVsOffsetData_HighResSpeed')
save('speedVsOffsetData_FieldsPerCycleSpeedUnits')
%%
zeroTimePhaseEdges=linspace(0,360,361);
zeroTimePhaseBins=edgesToBins(zeroTimePhaseEdges);

lowSpeedZeroTimeIdxes=abs(allLowSpeedPhaseVsTimeFracs(:,1)-0)<0.001;
midSpeedZeroTimeIdxes=abs(allMidSpeedPhaseVsTimeFracs(:,1)-0)<0.001;
highSpeedZeroTimeIdxes=abs(allHighSpeedPhaseVsTimeFracs(:,1)-0)<0.001;
highestSpeedZeroTimeIdxes=abs(allHighestSpeedPhaseVsTimeFracs(:,1)-0)<0.001;


lowSpeedZeroTimePhases=allLowSpeedPhaseVsTimeFracs(lowSpeedZeroTimeIdxes,2);
midSpeedZeroTimePhases=allMidSpeedPhaseVsTimeFracs(midSpeedZeroTimeIdxes,2);
highSpeedZeroTimePhases=allHighSpeedPhaseVsTimeFracs(highSpeedZeroTimeIdxes,2);
highestSpeedZeroTimePhases=allHighestSpeedPhaseVsTimeFracs(highestSpeedZeroTimeIdxes,2);


pl=getCircKernelDistr(lowSpeedZeroTimePhases,zeroTimePhaseEdges);
pm=getCircKernelDistr(midSpeedZeroTimePhases,zeroTimePhaseEdges);
ph=getCircKernelDistr(highSpeedZeroTimePhases,zeroTimePhaseEdges);
phst=getCircKernelDistr(highestSpeedZeroTimePhases,zeroTimePhaseEdges);

%{
Nl=histcounts(lowSpeedZeroTimePhases,zeroTimePhaseEdges);
Nm=histcounts(midSpeedZeroTimePhases,zeroTimePhaseEdges);
Nh=histcounts(highSpeedZeroTimePhases,zeroTimePhaseEdges);
Nhst=histcounts(highestSpeedZeroTimePhases,zeroTimePhaseEdges);
%}

figure; 
pL=plotCircXDistr(zeroTimePhaseBins,pl,'k-')
hold on
pM=plotCircXDistr(zeroTimePhaseBins,pm,'b-')
pH=plotCircXDistr(zeroTimePhaseBins,ph,'g-')
pHST=plotCircXDistr(zeroTimePhaseBins,phst,'r-')

%{
plot(zeroTimePhaseBins,Nnormalizedum(Nm),'b-','LineWidth',3)
plot(zeroTimePhaseBins,Nh/sum(Nh),'g-','LineWidth',3)
plot(zeroTimePhaseBins,Nhst/sum(Nhst),'r-','LineWidth',3)
%}

xlabel('Zero-time phase')
ylabel('Probability')


plot([360 360],ylim,'k--','LineWidth',2)
plot([0 0],ylim,'k--','LineWidth',2)
plot([720 720],ylim,'k--','LineWidth',2)
legend([pL,pM,pH,pHST],{sprintf('speed < %.2f normalized',singleCellLowSpeedMax), sprintf('%.2f normalized < speed < %.2f normalized',singleCellLowSpeedMax,singleCellMidSpeedMin),sprintf('%.2f normalized < speed < %.2f normalized',singleCellMidSpeedMin,singleCellHighSpeedMin),...
    sprintf('speed > %.2f normalized',singleCellHighSpeedMin)})

legend boxoff
setFigFontTo(16)
%%
figure

%pPhaseDiffEdges=linspace(-180,180,31);
pPhaseDiffEdges=linspace(-360,360,41);
pPhaseDiffBinCenters=edgesToBins(pPhaseDiffEdges);

%[p] = getProbDist(allHighVsLow,pPhaseDiffEdges,1,1);
[p] = getProbDist(allHighVsLow,pPhaseDiffEdges,1,0);
figure;
plot(pPhaseDiffBinCenters,p,'r-','LineWidth',3)
hold on; plot([0 0], ylim,'k--','LineWidth',3)
box off

xlabel('High speed minus low speed avg phase difference (degrees)')
ylabel('Probability')
%xlim([-180 180])
xlim([-360 360])
ylim([-Inf Inf])

%meanPhaseDiff=circMeanDeg(allHighVsLow);
%meanPhaseSEM=getSEMacrossRows(allHighVsLow);

%{
[h mu ul ll] = circ_mtest(deg2rad(allHighVsLow), 0);


muDeg=round(rad2deg(mu));
ulDeg=round(rad2deg(ul));
llDeg=round(rad2deg(ll));
%}

%title({sprintf('Distribution of speed induced phase shift across fields (n=%d)',length(allHighVsLow)),'significantly less than 0 by circ m test',sprintf('mean=%d%c, lower limit=%d%c, upper limit=%d%c, 95%% confidence',muDeg,char(176),llDeg,char(176),ulDeg,char(176))})
setFigFontTo(18)
maxFigHalfWidth
saveas(gcf,'distrPhaseShiftBySpeedAcrossCells.png')
%%
phaseAdvanceFig=figure
phaseDelayFig=figure
plotPolar=0;
if(plotPolar)
polarPlot=figure;
end
cartesianPlot=figure
speedGroups=[2 1];
for ci=1:size(allCellPairedLowHigh,1)
    %if(abs(diff(allCellPairedLowHigh(ci,:)))<=180)% && allCellPairedLowHigh(ci,1) >180)
   if(angdiffDeg([allCellPairedLowHigh(ci,1) allCellPairedLowHigh(ci,2)]) < 0) 
        %if(abs(angdiffDeg([allCellPairedLowHigh(ci,1) allCellPairedLowHigh(ci,2)])) < 90)
            %figure(cartesianPlot)
             figure(phaseAdvanceFig)
            if(abs(diff(allCellPairedLowHigh(ci,:)))<=180)
                plot(1:2,allCellPairedLowHigh(ci,:),'r-','LineWidth',0.75)
                hold on
                 plot(1:2,allCellPairedLowHigh(ci,:),'k.','MarkerSize',10)
            else
                plot(1:2,[allCellPairedLowHigh(ci,1) allCellPairedLowHigh(ci,2)-360] ,'r-','LineWidth',0.75)
                hold on
                 plot(1:2,[allCellPairedLowHigh(ci,1) allCellPairedLowHigh(ci,2)-360] ,'k.','MarkerSize',10)
            end
            
            if(plotPolar)
                figure(polarPlot)

                polarplot(deg2rad(allCellPairedLowHigh(ci,:)),speedGroups,'r-')
                hold on
                polarplot(deg2rad(allCellPairedLowHigh(ci,:)),speedGroups,'k.','MarkerSize',10)
            end
          else
             
            if(plotPolar)
                figure(polarPlot)
                 polarplot(deg2rad(allCellPairedLowHigh(ci,:)),speedGroups,'b-')
                 hold on
                 polarplot(deg2rad(allCellPairedLowHigh(ci,:)),speedGroups,'k.','MarkerSize',10)
            end
             
              figure(phaseDelayFig)
              if(abs(diff(allCellPairedLowHigh(ci,:)))<=180)
                plot(1:2,allCellPairedLowHigh(ci,:),'b-','LineWidth',0.75)
                hold on
                 plot(1:2,allCellPairedLowHigh(ci,:),'k.','MarkerSize',10)
              else
                   plot(1:2,[allCellPairedLowHigh(ci,1) allCellPairedLowHigh(ci,2)+360] ,'b-','LineWidth',0.75)
                hold on
                 plot(1:2,[allCellPairedLowHigh(ci,1) allCellPairedLowHigh(ci,2)+360] ,'k.','MarkerSize',10)
              end
        end
    %end
   %end
    hold on
end


fracAdvancedWithSpeed=sum(allHighVsLow<0)/length(allHighVsLow)

 figure(phaseAdvanceFig)
xlim([0.5 2.5])
ylim([0 360])
%ylim([-180 540])
xlabel('Running speed group')
ylabel('Average theta phase')
%title(sprintf('Variation within single fields across laps, n=%d fields',length(allHighVsLow)))
    title({sprintf('Variation within single fields across laps from speed < %.2f normalized to speed > %.2f normalized',singleCellLowSpeedMax,singleCellHighSpeedMin),sprintf('phase advance in %d out of %d fields',sum(allHighVsLow<0),length(allHighVsLow))})


setFigFontTo(18)
maxFig

saveas(gcf,'singleCellVariationWithSpeedAcrossLapsAdvance.png')

figure(phaseDelayFig)
xlim([0.5 2.5])
ylim([0 360])
%ylim([-180 540])
xlabel('Running speed group')
ylabel('Average theta phase')
%title(sprintf('Variation within single fields across laps, n=%d fields',length(allHighVsLow)))
    title({sprintf('Variation within single fields across laps from speed < %.2f cnormalized to speed > %.2f cnormalized',singleCellLowSpeedMax,singleCellHighSpeedMin),sprintf('phase delay in %d out of %d fields',sum(allHighVsLow>0),length(allHighVsLow))})


setFigFontTo(18)
maxFig

saveas(gcf,'singleCellVariationWithSpeedAcrossLapsDelay.png')


nanmean(allHighVsLow)
getSEMacrossRows(allHighVsLow)


%figure; histogram(allRangesPerTraversal)

if(plotPolar)
    figure(polarPlot)
    labelPolarPlot('Avg theta phase','Running speed group',2)

    title(sprintf('Variation within single fields across laps, phase advance in %d out of %d fields',sum(allHighVsLow<0),length(allHighVsLow)))

    setFigFontTo(18)
    maxFig

    saveas(gcf,'singleCellVariationWithSpeedAcrossLapsCirc.tif')
end



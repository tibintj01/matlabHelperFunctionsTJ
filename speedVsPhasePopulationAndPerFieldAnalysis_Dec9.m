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

speedCategoryBounds=0.35:0.025:0.65; %high res (2.5cm/s bin width )
speedCategoryBounds=0.35:0.02:0.65; %high res (2cm/s bin width )
speedCategoryBounds=0.4:0.02:0.7; %high res (2cm/s bin width )
speedCategoryBounds=0.35:0.01:0.7; %high res (2cm/s bin width )
speedCategoryBounds=0.4:0.01:0.65; %high res (2cm/s bin width )
speedCategoryBounds=0.3:0.02:0.65; %high res (2cm/s bin width )
%speedCategoryBounds=[0 0.25 0.35 0.45 0.55 0.65 0.75];
speedCategoryBounds=0.35:0.01:0.7; %high res (2cm/s bin width )
speedCategoryBounds=0.4:0.01:0.65; %high res (2cm/s bin width )
speedCategoryBounds=0.4:0.01:0.6; %high res (2cm/s bin width )
speedCategoryBounds=0.3:0.01:0.6; %high res (2cm/s bin width )
speedCategoryBounds=0.35:0.025:0.7; %high res (2cm/s bin width )
speedCategoryBounds=0.35:0.025:0.7; %high res (2cm/s bin width )
%speedCategoryBounds=0.35:0.05:0.7; %high res (2cm/s bin width )
speedCategoryBounds=0.35:0.05:0.7; %high res (2cm/s bin width )
%speedCategoryBounds=0.35:0.025:0.7; %high res (2cm/s bin width )
%speedCategoryBounds=0.2:0.15:0.8;

%speedCategoryBounds=0.30:0.03:0.7;
speedCategoryBounds=0.35:0.05:0.7;
speedCategoryBounds=0.35:0.05:0.75;
speedCategoryBounds=0.30:0.05:0.75;
speedCategoryBounds=0.25:0.1:0.75;
speedCategoryBounds=0.3:0.1:0.8;
%speedCategoryBounds=0.2:0.1:0.8;


%speedCategoryBounds=[0.2000    0.4000    0.5000    0.6000    0.7000    0.8000];


speedCategoryBounds=[0.3000    0.4000    0.5000    0.6000    0.7000    0.8000];

speedCategoryBounds=0.35:0.05:0.7;
speedCategoryBounds=0.2:0.075:0.8;

speedCategoryBounds=0.35:0.025:0.75;
speedCategoryBounds=0.35:0.05:0.75;

speedCategoryBounds=0.35:0.025:0.75;

%speedCategoryBounds=0.1:0.15:0.7; %4 categories, 5 bounds

restrictToSameSizeFields=1;


minFieldSizeThetaCycles=7; %not great
%minFieldSizeThetaCycles=8; %not bad
minFieldSizeThetaCycles=9;
maxFieldSizeThetaCycles=Inf;
   

speedCategoryBounds=[0.1 0.35 0.6 0.85];

numSpeedGroups=length(speedCategoryBounds)-1;


allPhaseVsTimePerSpeedGroupAllFields=cell(numSpeedGroups,1);
zeroTimePhasePerSpeedGpAllFields=cell(numSpeedGroups,1);


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

minRunSpeedDisp=minRunningSpeed;%m/s
maxRunSpeedDisp=maxRunningSpeed; %m/s

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
%maxSpeedRange=0.03;

%maxSpeedRange=0.05;
%maxSpeedRange=0.1;


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
                
                currFieldTimeBound=(currUnitStruct.cycleTimeFieldBoundPerDirPerField{fii,di}-1); %-1 becaues manually entered before changing first theta cycle from 1 to 0

                
                if(std(meanSpeedPerTraversal)<minSpeedStdvAcrossLaps)
                    continue
                end
                
                if(restrictToSameSizeFields)
                    if(currFieldTimeBound>maxFieldSizeThetaCycles || currFieldTimeBound<minFieldSizeThetaCycles)
                        continue
                    end
                end
                
                semSpeedPerTraversal=currFieldThetaData.allLapInFieldSpeedRangePerTraversal;
                
                allRangesPerTraversal=[allRangesPerTraversal;semSpeedPerTraversal(:)];
                %meanRangePerTraversal=
                
                 uniformSpeedPerTraversal=semSpeedPerTraversal<=maxSpeedRange;
                 
   
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 %TAKE ONLY LAPS WITH ~UNIFORM SPEED THROUGH FIELD
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 meanSpeedPerTraversal=meanSpeedPerTraversal(uniformSpeedPerTraversal);
                 meanPhasePerTraversal=meanPhasePerTraversal(uniformSpeedPerTraversal);
                 
                 if(isempty(meanSpeedPerTraversal) || isempty(meanPhasePerTraversal))
                     continue
                 end
                 
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
                 zeroTimePhasePerSpeedGpCurrField=NaN(numSpeedGroups,1);
                    currSpeedGpIdxes=cell(numSpeedGroups,1);
                    
                    %zeroPhaseMaxTime=0.15;
                    zeroPhaseMaxTime=0.2;
                     %zeroPhaseMaxTime=0.3;
                    minNumZeroPhases=3;
                    
                  for si=1:numSpeedGroups
                      currSpeedBinMin=speedCategoryBounds(si);
                      currSpeedBinMax=speedCategoryBounds(si+1);
                      
                      currSpeedGpIdxes{si}=meanSpeedPerTraversal>=currSpeedBinMin & meanSpeedPerTraversal<currSpeedBinMax;
                      
                      currSpeedGpPhaseVsTime=vertcat(currFieldThetaData.allLapInFieldAllTimeVsPhasePerTraversal{currSpeedGpIdxes{si}});
                      allPhaseVsTimePerSpeedGroupCurrField{si}=currSpeedGpPhaseVsTime;
                      
                      if(isempty(currSpeedGpPhaseVsTime))
                          continue
                      end
                      
                      zeroPhaseIdxesCurrGp=currSpeedGpPhaseVsTime(:,1)<zeroPhaseMaxTime;
                      
                      if(sum(zeroPhaseIdxesCurrGp)>=minNumZeroPhases)
                        zeroTimePhasePerSpeedGpCurrField(si)=circMeanDeg( currSpeedGpPhaseVsTime(zeroPhaseIdxesCurrGp,2));
                      end
                  end
            
              

                currFieldTimeBound=currUnitStruct.cycleTimeFieldBoundPerDirPerField{fii,di}-1; %-1 becaues manually entered before changing first theta cycle from 1 to 0
                
                currFieldTimeBound=currFieldTimeBound*0.9*0.9*0.95;
                
                allSpikesCurrTimeFieldFracs=currFieldThetaData.allLapSpikeCycleIDsFromFieldEntry/currFieldTimeBound;
                        
                goodLapSpikesCurrTimeFieldFracs=allSpikesCurrTimeFieldFracs(allSpikesGoodLapIdxes);
                goodLapSpikesCurrSpeeds=allSpikesCurrSpeeds(allSpikesGoodLapIdxes);
                goodLapSpikesCurrTimeFieldPhases=allSpikesCurrTimeFieldPhases(allSpikesGoodLapIdxes);
                
                allSpikeExpTimes=allSpikeExpTimes(allSpikesGoodLapIdxes);
              
                allPhases=goodLapSpikesCurrTimeFieldPhases(:);
               
              %figure; plot(allPhaseVsTimePerSpeedGroupCurrField{3}(:,1),allPhaseVsTimePerSpeedGroupCurrField{3}(:,2),'ko')
                
              
                %disp('')
                

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
                    xlabel('Mean speed in traversal (m/s)')
                    
                    
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
                            
                            if(~isnan(zeroTimePhasePerSpeedGpCurrField(si)))
                             zeroTimePhasePerSpeedGpAllFields{si}=[zeroTimePhasePerSpeedGpAllFields{si}; zeroTimePhasePerSpeedGpCurrField(si)];
                            end
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
                    %ylabel('Speed (m/s)')
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
                        xlabel('Speed (m/s)')
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
                    %ylabel('Speed (m/s)')
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
ylabel(cb,'Speed (m/s)')
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
xlabel('Avg running speed in field traversal (m/s)')

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
ylabel(cb,'Avg running speed in field traversal (m/s)')

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
xlabel('Avg running speed in field traversal (m/s)')
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


xlabel('Avg running speed in field traversal (m/s)')
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
%timeEdges=linspace(0,1,35);
%timeEdges=linspace(0,1,25);
phaseEdges=linspace(0,360,31);
%phaseEdges=linspace(0,360,21);
allPhasesHiVsLowSpeedFig=figure;
withSmooth=1;

%for si=1:numSpeedGroups
for si=numSpeedGroups:-1:1
  currSpeedBinMin=speedCategoryBounds(si);
  currSpeedBinMax=speedCategoryBounds(si+1);
                      

  currSpeedGroupPhaseVsTime=allPhaseVsTimePerSpeedGroupAllFields{si};
  
    subplot(1,numSpeedGroups,si)
    getJointDistr(currSpeedGroupPhaseVsTime(:,1),currSpeedGroupPhaseVsTime(:,2),timeEdges,phaseEdges,allPhasesHiVsLowSpeedFig,withSmooth);
    %title(sprintf('speed < %.2f m/s',singleCellLowSpeedMax))
    title(sprintf('%.2f m/s < speed < %.2f m/s',currSpeedBinMin,currSpeedBinMax))
    
    if(si==numSpeedGroups)
        climits=caxis;
    else
        %caxis([climits])
    end
    
    caxis([0 2.4e-3])

end

%{
subplot(1,4,2)
getJointDistr(allMidSpeedPhaseVsTimeFracs(:,1),allMidSpeedPhaseVsTimeFracs(:,2),timeEdges,phaseEdges,allPhasesHiVsLowSpeedFig);
title(sprintf('%.2f m/s < speed < %.2f m/s',singleCellLowSpeedMax,singleCellMidSpeedMin))

subplot(1,4,3)
getJointDistr(allHighSpeedPhaseVsTimeFracs(:,1),allHighSpeedPhaseVsTimeFracs(:,2),timeEdges,phaseEdges,allPhasesHiVsLowSpeedFig);
title(sprintf('%.2f m/s < speed < %.2f m/s',singleCellMidSpeedMin,singleCellHighSpeedMin))

subplot(1,4,4)
getJointDistr(allHighestSpeedPhaseVsTimeFracs(:,1),allHighestSpeedPhaseVsTimeFracs(:,2),timeEdges,phaseEdges,allPhasesHiVsLowSpeedFig);
title(sprintf('speed > %.2f m/s',singleCellHighSpeedMin))
%}
setFigFontTo(16)
%%
figure
%timeEdges=linspace(0,1,12);
timeEdges=linspace(0,1,7);
timeEdges=linspace(0,1,8);
%timeEdges=linspace(0,1,12);
%timeEdges=linspace(0,1,6);

speedGroupColors=getBiColorMap(numSpeedGroups);

minCircShiftPhase=20;

for si=1:numSpeedGroups
  currSpeedBinMin=speedCategoryBounds(si);
  currSpeedBinMax=speedCategoryBounds(si+1);
                      

  currSpeedGroupPhaseVsTime=allPhaseVsTimePerSpeedGroupAllFields{si};
  
  [timeBinCenters,currSpeedGpPhaseAvgValues] = getBinnedCircAverages(currSpeedGroupPhaseVsTime(:,1),currSpeedGroupPhaseVsTime(:,2),timeEdges);
  
  currPhases=currSpeedGpPhaseAvgValues;
  
  shiftIdxes=currPhases<minCircShiftPhase;
  
  currPhases(shiftIdxes)=currPhases(shiftIdxes)+360;
    plot(timeBinCenters,currPhases,'-','Color', speedGroupColors(si,:),'LineWidth',3)
    hold on
     plot(timeBinCenters,currPhases,'.','Color', speedGroupColors(si,:),'MarkerSize',40)

end
colormap(getBiColorMap)
cb=colorbar('north')

ylabel(cb,'running speed (m/s)')
caxis([min(speedCategoryBounds) max(speedCategoryBounds)])
 xlabel('Time in field (frac)')
 ylabel('Mean theta phase (degrees)')
 box off
 setFigFontTo(16)
 
 ylim([50 380])
  %ylim([60 360])

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


%legend([pL,pM,pH,pHST],{sprintf('speed < %.2f m/s',singleCellLowSpeedMax), sprintf('%.2f m/s < speed < %.2f m/s',singleCellLowSpeedMax,singleCellMidSpeedMin),sprintf('%.2f m/s < speed < %.2f m/s',singleCellMidSpeedMin,singleCellHighSpeedMin),sprintf('speed > %.2f m/s',singleCellHighSpeedMin)})
%legend([pL,pM,pH,pHST],{sprintf('speed < %.2f m/s',singleCellLowSpeedMax), sprintf('%.2f m/s < speed < %.2f m/s',singleCellLowSpeedMax,singleCellMidSpeedMin),sprintf('%.2f m/s < speed < %.2f m/s',singleCellMidSpeedMin,singleCellHighSpeedMin),sprintf('speed > %.2f m/s',singleCellHighSpeedMin)})

%legend boxoff

%legend([p1 p2 p3 p4], {sprintf('speed < %.2f m/s',singleCellLowSpeedMax),...
%    sprintf('%.2f m/s < speed < %.2f m/s',singleCellLowSpeedMax,singleCellMidSpeedMin),sprintf('%.2f m/s < speed < %.2f m/s',singleCellMidSpeedMin,singleCellHighSpeedMin), sprintf('speed > %.2f m/s',singleCellHighSpeedMin)})



%save('speedVsOffsetData_Oct8_2021.mat')
%save('speedVsOffsetData_CloseTimeBound')
%save('speedVsOffsetData_HighResSpeed_Oct25')
save('speedVsOffsetData_HighResSpeed_Nov5')
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
plot(zeroTimePhaseBins,Nm/sum(Nm),'b-','LineWidth',3)
plot(zeroTimePhaseBins,Nh/sum(Nh),'g-','LineWidth',3)
plot(zeroTimePhaseBins,Nhst/sum(Nhst),'r-','LineWidth',3)
%}

xlabel('Zero-time phase')
ylabel('Probability')


plot([360 360],ylim,'k--','LineWidth',2)
plot([0 0],ylim,'k--','LineWidth',2)
plot([720 720],ylim,'k--','LineWidth',2)
legend([pL,pM,pH,pHST],{sprintf('speed < %.2f m/s',singleCellLowSpeedMax), sprintf('%.2f m/s < speed < %.2f m/s',singleCellLowSpeedMax,singleCellMidSpeedMin),sprintf('%.2f m/s < speed < %.2f m/s',singleCellMidSpeedMin,singleCellHighSpeedMin),...
    sprintf('speed > %.2f m/s',singleCellHighSpeedMin)})

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
    title({sprintf('Variation within single fields across laps from speed < %.2f m/s to speed > %.2f m/s',singleCellLowSpeedMax,singleCellHighSpeedMin),sprintf('phase advance in %d out of %d fields',sum(allHighVsLow<0),length(allHighVsLow))})


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
    title({sprintf('Variation within single fields across laps from speed < %.2f cm/s to speed > %.2f cm/s',singleCellLowSpeedMax,singleCellHighSpeedMin),sprintf('phase delay in %d out of %d fields',sum(allHighVsLow>0),length(allHighVsLow))})


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



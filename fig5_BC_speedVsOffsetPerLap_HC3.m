 
close all; clear all; clc
%dataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles';
processedDataDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData';
unitDataDir=fullfile(processedDataDir,'unitSpikeInfoStructs');
dataDir=unitDataDir;

filePaths=getFilePathsRegex(dataDir,'*mat');


masterSpeeds=[];
masterCircMeanPhase=[];
masterMeanTimeIn=[];
masterMeanTimeTo=[];
masterTimeSinceExpStart=[];

masterSpeedMeanPhaseCircR=[];
masterNumLaps=[];
masterSpeedVar=[];
masterFieldIDs=[];
masterBinnedPhaseAvgs=[];
masterSpeedBins=[];

numSpeedBins=20;
numSpeedBins=30;
%numSpeedBins=10;
speedEdges=linspace(0,1.2,numSpeedBins+1);
speedBins=edgesToBins(speedEdges);
%numSpeedBins=length(speedBins);

showEachFieldPlots=1;
showEachFieldPlots=0;


avgPerSpeedBeforeCombiningPop=1;
avgPerSpeedBeforeCombiningPop=0;

minNumLaps=8;
threshSpeedVarAcrossLaps=0.3;
threshSpeedVarAcrossLaps=0.2;
threshSpeedVarAcrossLaps=0.1;
threshSpeedVarAcrossLaps=0.2;
threshSpeedVarAcrossLaps=0.15;
threshSpeedVarAcrossLaps=0;
threshSpeedVarAcrossLaps=0.2;
threshSpeedVarAcrossLaps=0.25;
threshSpeedVarAcrossLaps=0.18;
%threshSpeedVarAcrossLaps=0.2; 
threshSpeedVarAcrossLaps=0.2328; %prctile(masterSpeeds,25)
minSpeedVar=0.15;
%minSpeedVar=0.1; %dip in speed var histogram
minSpeedVar=0; %dip in speed var histogram

%{

speedRangeFieldData=load('fieldsWithSpeedRange.mat');

fieldIDsWithSpeedRangeAcrossLaps=speedRangeFieldData.fieldIDsWithLargeSpeedRange;
%}

fieldCount=0;
dispFieldCount=0;
for i=1:length(filePaths)
    i
    currFilePath=filePaths{i};
    data=load(currFilePath);
    

    for di=1:2
        if(di==1)
            currDirStr='rightward';
            currStatsFieldName=sprintf('speedVsMeanPhaseStatsLeftwardNoAutoShiftRef');
        else
               currDirStr='leftward';
               currStatsFieldName=sprintf('speedVsMeanPhaseStatsLeftwardNoAutoShiftRef');
        end
        
        if(~isfield(data,(currStatsFieldName)))
            continue
        end
        
        if(~isfield(data.(currStatsFieldName),currDirStr))
            continue
        end
   
        
        numFields=data.(currStatsFieldName).(currDirStr).currNumFields;
        
        for fii=1:numFields
        %for fii=1:1
           
            %avgFieldSpeedPerLap= data.(currStatsFieldName).(currDirStr).avgFieldSpeedPerLap(fii,:);
            avgFieldSpeedPerLap=abs(data.(currStatsFieldName).(currDirStr).avgFieldSpeedPerLap(fii,:));
            
            
            circAvgPhaseInFieldPerLap= data.(currStatsFieldName).(currDirStr).circAvgPhaseInFieldPerLap(fii,:);
            
            avgFieldSpeedPerLap=avgFieldSpeedPerLap(1:length(circAvgPhaseInFieldPerLap));
            
            
            if(nanstd(avgFieldSpeedPerLap)<threshSpeedVarAcrossLaps)
                continue
            end
            
            
       
             %tooSlowIdxes=avgFieldSpeedPerLap<0.05;
            
            %badIdxes=avgFieldSpeedPerLap<0.05 | avgFieldSpeedPerLap> 1.2;
            badIdxes=avgFieldSpeedPerLap<0 | avgFieldSpeedPerLap> 1.2;
            %badIdxes=avgFieldSpeedPerLap<0.5 | avgFieldSpeedPerLap> 1.2;
            
            avgTimeInFieldPerLap= data.(currStatsFieldName).(currDirStr).avgTimeInFieldPerLap(fii,:);
            avgTimeToFieldEndPerLap= data.(currStatsFieldName).(currDirStr).avgTimeToFieldEndPerLap(fii,:);
            
            
            %badIdxes=badIdxes | circAvgPhaseInFieldPerLap<180;
            circMedianPhaseInFieldPerLap= data.(currStatsFieldName).(currDirStr).circMedianPhaseInFieldPerLap(fii,:);
            
            timeSinceExpStartPerLap=data.(currStatsFieldName).(currDirStr).timeSinceExpStartPerLap(fii,:);
            
            %circAvgPhaseInFieldPerLap=circMedianPhaseInFieldPerLap;
            %phaseTimeLinFitOffsetInFieldPerLap= data.(currStatsFieldName).(currDirStr).phaseTimeLinFitOffsetInFieldPerLap(fii,:);
            %phaseTimeLinFitSlopeInFieldPerLap=data.(currStatsFieldName).(currDirStr).phaseTimeLinFitSlopeInFieldPerLap(fii,:);

            circAvgPhaseInFieldPerLap(badIdxes)=[];
            avgTimeInFieldPerLap(badIdxes)=[];
            avgTimeToFieldEndPerLap(badIdxes)=[];
            
            circMedianPhaseInFieldPerLap(badIdxes)=[];
            %phaseTimeLinFitOffsetInFieldPerLap(badIdxes)=[];
            %phaseTimeLinFitSlopeInFieldPerLap(badIdxes)=[];
            avgFieldSpeedPerLap(badIdxes)=[];
            
            timeSinceExpStartPerLap(badIdxes)=[];
            
            
                 
                nonNanIdxes=~(isnan(circAvgPhaseInFieldPerLap) | isnan(avgFieldSpeedPerLap));
                 numNonNanPts=sum(nonNanIdxes);
                 %numLaps=length(avgFieldSpeedPerLap);
                 numLaps=numNonNanPts;
                 if(numNonNanPts<minNumLaps)
                     continue
                 end
                 timeSinceExpStartPerLap=timeSinceExpStartPerLap(nonNanIdxes);
                 circAvgPhaseInFieldPerLap=circAvgPhaseInFieldPerLap(nonNanIdxes);
                 avgFieldSpeedPerLap=avgFieldSpeedPerLap(nonNanIdxes);
                 avgTimeInFieldPerLap=avgTimeInFieldPerLap(nonNanIdxes);
                 avgTimeToFieldEndPerLap=avgTimeToFieldEndPerLap(nonNanIdxes);
                 
                 %{
                  currSpeedVar=nanstd(avgFieldSpeedPerLap(:));
            if(currSpeedVar<threshSpeedVarAcrossLaps)
                continue
            end
                 %}
                 
                 
                 %{
                     numShifts=360;
                     varArray=NaN(numShifts,1);
                     for shiftIdx=1:numShifts

                         currSampleMeanPhasesShifted=mod(circAvgPhaseInFieldPerLap+shiftIdx,360);

                         varArray(shiftIdx)=nanstd(currSampleMeanPhasesShifted);
                     end
                     [~,bestShift]=min(varArray);
                     circAvgPhaseInFieldPerLap=mod(circAvgPhaseInFieldPerLap+bestShift,360);
                 %}
            
            commonXLims=[0 1.2];
            commonYLims=[0 360];
               
            fieldCount=fieldCount+1
            
            %{
            if(min(abs(fieldIDsWithSpeedRangeAcrossLaps-fieldCount))>0)
                %continue
            end
            %}
            
            circAvgPhaseInFieldPerLap=circAvgPhaseInFieldPerLap(:);
            avgFieldSpeedPerLap=avgFieldSpeedPerLap(:);
            notNanIdx=~isnan(avgFieldSpeedPerLap(:))&(~isnan(circAvgPhaseInFieldPerLap(:)));
            
            %{
            [sortedSpeedPerLap,sortIdx]=sort(avgFieldSpeedPerLap);
            sortedCircAvgPhaseInFieldPerLap=circAvgPhaseInFieldPerLap(sortIdx);
            
            ybar = scatstat1(sortedSpeedPerLap,sortedCircAvgPhaseInFieldPerLap,0.1);
            %}
            discretizedSpeeds=discretize(avgFieldSpeedPerLap,speedEdges);

                avgPhaseInSpeedBin=NaN(numSpeedBins,1);
                for bi=1:numSpeedBins
                    currSpeedBinIdxes=discretizedSpeeds==bi;

                    if(sum(currSpeedBinIdxes)>1)
                    avgPhaseInSpeedBin(bi)=circMeanDeg(circAvgPhaseInFieldPerLap(currSpeedBinIdxes));
                    end

                end
            
             if(showEachFieldPlots)
                 dispFieldCount=dispFieldCount+1;
                
                  figure(1);
                %subplot(14,13,fieldCount)
                subplot(11,11,fieldCount)
                %subplot(6,6,dispFieldCount)
                
                
                %plot(avgFieldSpeedPerLap,circAvgPhaseInFieldPerLap,'k.')
                %{
                plot(avgFieldSpeedPerLap,circAvgPhaseInFieldPerLap,'k.')
                hold on
                plot(speedBins,avgPhaseInSpeedBin,'ro')
                %}

                     ylim(commonYLims)
                     xlim(commonXLims)
                     %title(sprintf('std %.2f, %d',nanstd(avgFieldSpeedPerLap),fieldCount))
                     title(sprintf('std %.2f, %d',nanstd(avgFieldSpeedPerLap),fieldCount))
                     maxFig
                 
                     
            end
            
                 %{
             subplot(2,2,2)
            plot(avgFieldSpeedPerLap,circMedianPhaseInFieldPerLap,'ko')
                 ylim(commonYLims)
                  xlim(commonXLims)
             subplot(2,2,3)
            plot(avgFieldSpeedPerLap,phaseTimeLinFitOffsetInFieldPerLap,'ko')
                 %ylim(commonYLims)
             subplot(2,2,4)
            plot(avgFieldSpeedPerLap,phaseTimeLinFitSlopeInFieldPerLap,'ko')
            %}
            
            masterBinnedPhaseAvgs=[masterBinnedPhaseAvgs;avgPhaseInSpeedBin(:)];
            masterSpeedBins=[masterSpeedBins;speedBins(:)];
            
            masterSpeeds=[masterSpeeds; avgFieldSpeedPerLap(:)];
            
            masterCircMeanPhase=[masterCircMeanPhase; circAvgPhaseInFieldPerLap(:)];
            
             masterMeanTimeIn=[masterMeanTimeIn; avgTimeInFieldPerLap(:)];
              masterMeanTimeTo=[masterMeanTimeTo; avgTimeToFieldEndPerLap(:)];
              
              masterTimeSinceExpStart=[masterTimeSinceExpStart; timeSinceExpStartPerLap(:)];
              
            
            
            %fieldCount=fieldCount+1;
            masterFieldIDs=[masterFieldIDs; repelem(fieldCount,numLaps)'];
                
            
            masterNumLaps=[masterNumLaps; numLaps];
            %masterSpeedVar=[masterSpeedVar;currSpeedVar];
            
            %[ rho,p,s,b ] = kempter_lincirc( avgFieldSpeedPerLap(:),circAvgPhaseInFieldPerLap(:)*(pi/180) )
            %[ rho,p,s,b ] = kempter_lincirc( avgFieldSpeedPerLap(:),phaseTimeLinFitOffsetInFieldPerLap(:)*(pi/180) );
              %daspect([1 360 1])
            
            %corrParameter=masterCircMeanPhasel
            %{
            notNanIdx=~isnan(masterSpeeds)&(~isnan(masterCircMeanPhase(:)));
            rMat=corrcoef([masterSpeeds(notNanIdx),masterCircMeanPhase(notNanIdx)]);
            rho=rMat(1,2);
            %}
            
            
            
            %rMat=corrcoef([avgFieldSpeedPerLap(notNanIdx),circAvgPhaseInFieldPerLap(notNanIdx)]);
            %rho=rMat(1,2)
            
            
            %masterSpeedMeanPhaseCircR=[masterSpeedMeanPhaseCircR; rho];
            %ylim(commonYLims)
            %drawnow
        end
    
    end
end
%%

%masterHighSpeedIdxes=masterSpeeds>0.7;
%masterLowSpeedIdxes=masterSpeeds<0.2;

%masterHighSpeedIdxes=masterSpeeds>1;

%masterFieldIDs


minHighSpeed=0.6;
maxLowSpeed=0.2;

minHighSpeed=0.6;
maxLowSpeed=0.3;

minHighSpeed=0.5;
maxLowSpeed=0.5;
minHighSpeed=0.6;
maxLowSpeed=0.3;
minHighSpeed=0.6;
maxLowSpeed=0.2;

minHighSpeed=0.7;
maxLowSpeed=0.15;

minHighSpeed=0.6;
maxHighSpeed=1.2;

maxLowSpeed=0.3;
minLowSpeed=0.05;

minHighSpeed=0.8;
maxHighSpeed=1.2;

maxLowSpeed=0.2;
minLowSpeed=0.05;

minHighSpeed=0.7;
maxHighSpeed=1.2;

maxLowSpeed=0.35;
minLowSpeed=0.05;

minHighSpeed=0.6;
maxHighSpeed=1.2;

maxLowSpeed=0.3;
minLowSpeed=0.05;

%{
minHighSpeed=0.55;
maxHighSpeed=0.65;

maxLowSpeed=0.35;
minLowSpeed=0.25;
%}

if(avgPerSpeedBeforeCombiningPop)
    masterHighSpeedIdxes=masterSpeedBins>minHighSpeed;
    masterLowSpeedIdxes=masterSpeedBins<maxLowSpeed;
    masterCircMeanPhaseHighSpeeds=masterBinnedPhaseAvgs(masterHighSpeedIdxes)/360;
    masterCircMeanPhaseLowSpeeds=masterBinnedPhaseAvgs(masterLowSpeedIdxes)/360;
else
    masterHighSpeedIdxes=masterSpeeds>minHighSpeed & masterSpeeds<maxHighSpeed ;
    masterLowSpeedIdxes=masterSpeeds<maxLowSpeed & masterSpeeds>minLowSpeed ;

    masterCircMeanPhaseHighSpeeds=masterCircMeanPhase(masterHighSpeedIdxes)/360;
    masterCircMeanPhaseLowSpeeds=masterCircMeanPhase(masterLowSpeedIdxes)/360;
    
end

masterMeanTimeToHighSpeeds=masterMeanTimeTo(masterHighSpeedIdxes);
masterMeanTimeToLowSpeeds=masterMeanTimeTo(masterLowSpeedIdxes);
masterTimeSinceExpStartHighSpeeds=masterTimeSinceExpStart(masterHighSpeedIdxes);
masterTimeSinceExpStartLowSpeeds=masterTimeSinceExpStart(masterLowSpeedIdxes);

onlyLater=1;
%onlyLater=0;
if(onlyLater)
    minNoveltyThresh=5*60; %sec
     minNoveltyThresh=10*60; %sec
     minNoveltyThresh=15*60; %sec
     minNoveltyThresh=20*60; %sec
    masterPost10MinIdxesHigh=masterTimeSinceExpStartHighSpeeds>=minNoveltyThresh;
    masterPost10MinIdxesLow=masterTimeSinceExpStartLowSpeeds>=minNoveltyThresh;

    masterCircMeanPhaseHighSpeeds=masterCircMeanPhaseHighSpeeds(masterPost10MinIdxesHigh);
    masterCircMeanPhaseLowSpeeds=masterCircMeanPhaseLowSpeeds(masterPost10MinIdxesLow);
    masterTimeSinceExpStartHighSpeeds=masterTimeSinceExpStartHighSpeeds(masterPost10MinIdxesHigh);
    masterTimeSinceExpStartLowSpeeds=masterTimeSinceExpStartLowSpeeds(masterPost10MinIdxesLow);
end


%%
 figure; 
 numNoveltyBins=50;
 noveltyTimeEdges=linspace(0,2800,numNoveltyBins+1);
 
 %{
 [pNoveltyHighSpeed] = getProbDist(masterTimeSinceExpStartHighSpeeds,noveltyTimeEdges,1,0); %showPlots,isCirc
 hold on
 [pNoveltyLowSpeed] = getProbDist(masterTimeSinceExpStartLowSpeeds,noveltyTimeEdges,1,0);
 %}
 
  [pNoveltyHighSpeed] = getProbDist(masterTimeSinceExpStartHighSpeeds,noveltyTimeEdges,0,0); %showPlots,isCirc
 hold on
 [pNoveltyLowSpeed] = getProbDist(masterTimeSinceExpStartLowSpeeds,noveltyTimeEdges,0,0);
 disp('')
 
 
 % hold on; histogram(masterTimeSinceExpStartLowSpeeds,50)
 



%{
masterMeanTimeInHighSpeeds=masterMeanTimeIn(masterHighSpeedIdxes);
masterMeanTimeInLowSpeeds=masterMeanTimeIn(masterLowSpeedIdxes);


masterMeanTimeToHighSpeeds=masterMeanTimeTo(masterHighSpeedIdxes);
masterMeanTimeToLowSpeeds=masterMeanTimeTo(masterLowSpeedIdxes);
%}

figure;
%phaseEdges=0:0.015:1;
timeEdges=0:0.05:5;
timeEdges=0:0.025:5;
timeToBins=edgesToBins(timeEdges);

[Ntlow,c]=histcounts(masterMeanTimeToLowSpeeds,timeEdges);
hold on
[Nthigh,c]=histcounts(masterMeanTimeToHighSpeeds,timeEdges);

pTimeSlow=(Ntlow/sum(Ntlow(:)));
pTimeFast=(Nthigh/sum(Nthigh(:)));

smoothWind=15;
%smoothWind=10;
%smoothWind=5;
pTimeSlow=smooth(pTimeSlow,smoothWind);
pTimeFast=smooth(pTimeFast,smoothWind);

maxTime=5;
ptSlow=plot(maxTime-timeToBins,pTimeSlow,'b-','LineWidth',5)
hold on

ptFast=plot(maxTime-timeToBins,pTimeFast,'r-','LineWidth',5)

xlabel('Time to end of field (sec)')
ylabel('Probability')
title(sprintf('All place fields with speed std dev across laps > %.2f m/s',threshSpeedVarAcrossLaps))
legend([ptSlow ptFast],sprintf('running speed < %.2f m/s (n=%d samples, %d fields)',maxLowSpeed,sum(Ntlow(:)),length(masterNumLaps)),sprintf('running speed > %.2f m/s (n=%d samples,%d fields)',minHighSpeed,sum(Nthigh(:)),length(masterNumLaps)))


%xlim([0 360])
%xlim([60 300])
setFigFontTo(32)
maxFig




figure;
%phaseEdges=0:0.015:1;
phaseEdges=0:0.02:1;
phaseBinWidth=0.025;
phaseEdges=0:phaseBinWidth:1;

%phaseEdges=phaseEdges*360;

[Nlow,c]=histcounts(masterCircMeanPhaseLowSpeeds,phaseEdges);
hold on
[Nhigh,c]=histcounts(masterCircMeanPhaseHighSpeeds,phaseEdges);

pPhaseSlow=(Nlow/sum(Nlow(:)));
pPhaseFast=(Nhigh/sum(Nhigh(:)));

allFieldMeanPhasesLowSpeed=masterCircMeanPhaseLowSpeeds*360;
allFieldMeanPhasesHighSpeed=masterCircMeanPhaseHighSpeeds*360;
phaseEdgesDeg=phaseEdges*360;
phaseBinsDeg=edgesToBins(phaseEdgesDeg);
phaseBinWidthDeg=phaseBinWidth*360;
numPhaseBins=length(phaseBinsDeg);

save('lowSpeedVsHighSpeedPhaseProbDist.mat','masterMeanTimeToLowSpeeds', 'masterMeanTimeToHighSpeeds','pPhaseSlow','pPhaseFast','allFieldMeanPhasesLowSpeed','allFieldMeanPhasesHighSpeed','phaseEdgesDeg','phaseBinsDeg','phaseBinWidthDeg','numPhaseBins')


pPhaseSlow=[pPhaseSlow(:); pPhaseSlow(:); pPhaseSlow(:); pPhaseSlow(:)];
pPhaseFast=[pPhaseFast(:); pPhaseFast(:); pPhaseFast(:); pPhaseFast(:)];

smoothWind=15;
%smoothWind=10;
%smoothWind=5;
pPhaseSlow=smooth(pPhaseSlow,smoothWind);
pPhaseFast=smooth(pPhaseFast,smoothWind);
startBinPhase=numPhaseBins+1;
endBinPhase=startBinPhase+numPhaseBins-1;

pPhaseSlow=pPhaseSlow(startBinPhase:endBinPhase);
pPhaseFast=pPhaseFast(startBinPhase:endBinPhase);

pPhaseSlow=pPhaseSlow/sum(pPhaseSlow(:));
pPhaseFast=pPhaseFast/sum(pPhaseFast(:));

%degreesAxis=[edgesToBins(phaseEdges)*360-360, edgesToBins(phaseEdges)*360, edgesToBins(phaseEdges)*360+360 , edgesToBins(phaseEdges)*360+360*2];
degreesAxis=[edgesToBins(phaseEdges)*360];


pSlow=plot(degreesAxis,pPhaseSlow,'b-','LineWidth',5)
hold on

pFast=plot(degreesAxis,pPhaseFast,'r-','LineWidth',5)
xlim([-360 720])

plot([0 0],ylim,'k--','LineWidth',5)
plot([0 0]-360,ylim,'k--','LineWidth',5)
plot([0 0]+360,ylim,'k--','LineWidth',5)
plot([0 0]+360*2,ylim,'k--','LineWidth',5)

xlabel('Mean theta phase per lap (degrees)')
ylabel('Probability')
title(sprintf('All place fields with speed std dev across laps > %.2f m/s',threshSpeedVarAcrossLaps))
legend([pSlow pFast],sprintf('running speed < %.2f m/s (n=%d samples, %d fields)',...
    maxLowSpeed,sum(Nlow(:)),length(masterNumLaps)),sprintf('running speed > %.2f m/s (n=%d samples,%d fields)',...
    minHighSpeed,sum(Nhigh(:)),length(masterNumLaps)),'Location','best')

xlim([0 360])
%xlim([60 300])
ylim([0.007 0.047])
setFigFontTo(32)
maxFig

lowSpeedAnglesRad=2*pi*masterCircMeanPhaseLowSpeeds;
highSpeedAnglesRad=2*pi*masterCircMeanPhaseHighSpeeds;
%[pval table] = circ_wwtest(lowSpeedAnglesRad,highSpeedAnglesRad);

nanmean(masterTimeSinceExpStartHighSpeeds)
nanmean(masterTimeSinceExpStartLowSpeeds)

getSEMacrossRows(masterTimeSinceExpStartHighSpeeds)
getSEMacrossRows(masterTimeSinceExpStartLowSpeeds)

%[h,p,ci,stats] = ttest2(masterTimeSinceExpStartHighSpeeds,masterTimeSinceExpStartLowSpeeds)
[p,h,stats] = ranksum(masterTimeSinceExpStartHighSpeeds,masterTimeSinceExpStartLowSpeeds)

%%
runBootstrapTest=0;
if(runBootstrapTest)
    figure
    highSpeedVsLowSpeedPhaseKLDiv
end



%{
%%
fH=figure;
goodIdxes=masterSpeeds<1.2;
plotJointHeatMap(masterSpeeds(goodIdxes),masterCircMeanPhase(goodIdxes)/360,0.02,fH)
ylim([0.3 1])



%%
figure; 
histogram(masterSpeedVar)

bootstrapSampleSize=ceil(nanmean(masterNumLaps))
numShuffles=10000;

randomSampleRs=NaN(numShuffles,1);

for si=1:numShuffles
    if(mod(si,1000)==0)
        disp(si/numShuffles)
    end
    sampleSpeedIdxes=randsample(length(masterSpeeds),bootstrapSampleSize);
    currSampleSpeeds=masterSpeeds(sampleSpeedIdxes);
    
    sampleMeanPhaseIdxes=randsample(length(masterCircMeanPhase),bootstrapSampleSize);
    currSampleMeanPhases=masterCircMeanPhase(sampleMeanPhaseIdxes);
     %currSampleMeanPhases=masterCircMeanPhase(sampleSpeedIdxes);
    
     notNanIdx=~isnan(currSampleSpeeds)&(~isnan(currSampleMeanPhases(:)));
     %{
     numShifts=360;
     varArray=NaN(numShifts,1);
     for shiftIdx=1:numShifts
         
         currSampleMeanPhasesShifted=mod(currSampleMeanPhases+shiftIdx,360);
         
         varArray(shiftIdx)=nanstd(currSampleMeanPhasesShifted);
     end
     [~,bestShift]=min(varArray);
     currSampleMeanPhases=mod(currSampleMeanPhases+bestShift,360);
     %}
     
            rMat=corrcoef([currSampleSpeeds(notNanIdx),currSampleMeanPhases(notNanIdx)]);
            rho=rMat(1,2);
            
            randomSampleRs(si)=rho;
end
%%

figure; 
speedVsMeanPhaseInFieldEdges=-0.1:0.005:0.1;
speedVsMeanPhaseInFieldEdges=-0.5:0.02:0.5;
[Nreal,c]=histcounts(masterSpeedMeanPhaseCircR,speedVsMeanPhaseInFieldEdges);

hold on
[Nboot,c]=histcounts(randomSampleRs,speedVsMeanPhaseInFieldEdges);
pReal=Nreal/sum(Nreal);
pBoot=Nboot/sum(Nboot);

binCenters=edgesToBins(speedVsMeanPhaseInFieldEdges);

plot(binCenters,pReal,'r-','LineWidth',5)
hold on
plot(binCenters,pBoot,'k-','LineWidth',5)
%{
%%
figure; plot(masterSpeeds,masterCircMeanPhase,'k.')
%%
[ rho,p,s,b ] = kempter_lincirc( masterSpeeds,masterCircMeanPhase*(pi/180) )
%title(sprintf('circ-linear R=%.3f',rho))

binWidth=0.01;
fH=figure
inBoundIdx=masterSpeeds<1.2;
plotJointHeatMap(masterSpeeds(inBoundIdx),masterCircMeanPhase(inBoundIdx)/360,binWidth,fH)

%}

%}
%}

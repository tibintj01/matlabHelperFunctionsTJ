%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load unit struct file names
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; clc
processedDataDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData';
unitDataDir=fullfile(processedDataDir,'unitSpikeInfoStructs');
exampleFieldsDir='exampleFieldsPhasePrecessionLowVsHighSpeed_ZtoMinus3_ConstantSpeedTravs';

touchDir(exampleFieldsDir)

saveLFPcycleInfoDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData/lfpCycleInfos';
saveRawLFPInfoDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData/rawLFPstructs';

filePaths=getFilePathsRegex(unitDataDir,'*mat');

spikeDataPerField=load('allUnitInFieldGoodSpeedSpikePhasesAndTimes.mat');

recalcSpeedperSpike=0;

recalcSpeedperSpike=1;
recalcThetaData=1;


recalcSpeedperSpike=0;
recalcThetaData=0;

justExamplePlots=1;


containsSpikeTol=0.001;

maxSpeedStdDev=0.1;
%maxSpeedStdDev=0.05;

minDataPerBin=2;

showPlots=0;

excludeHighVarTraversals=0;

if(recalcSpeedperSpike)

     %spikeSpeedsPerField=cell(totalFieldCount,1);

    spikeTimesInExpPerField=spikeDataPerField.spikeTimesInExpPerField;
    spikePhasesPerField =spikeDataPerField.spikePhasesPerField;
    unitInfoPathPerField=spikeDataPerField.unitInfoPathPerField;
    spikeSpeedsPerField=spikeDataPerField.spikeSpeedsPerField;
    dirPerField=spikeDataPerField.dirPerField;
    
    totalFieldCount=spikeDataPerField.totalFieldCount;
    
    if(~recalcThetaData)
        thetaDataPerField=spikeDataPerField.thetaDataPerField;
    else
          thetaDataPerField=cell(totalFieldCount,1);
          fieldNumWithinUnitPerField=cell(totalFieldCount,1);
    end
    
        
        
    
    spikeAvgSpeedInTraversalsPerField=cell(totalFieldCount,1);


    startFieldIdx=1;
     fH=figure;
        
   
        
    for fi=startFieldIdx:totalFieldCount
        %for fi=76:length(filePaths)
        tic
           fi
        currFilePath=unitInfoPathPerField{fi};
        currFileName=getFileNameFromPath(currFilePath);
        fileBaseName=currFileName(1:(end-4));

        dataStruct=load(currFilePath);
        data=dataStruct.unitInfo;
        
        currFieldDirNum=dirPerField{fi};
        
        if(currFieldDirNum==1)
            currFieldDirStr='rightward';
        elseif(currFieldDirNum==2)
            currFieldDirStr='leftward';
        end
        
        if(~recalcThetaData)
             currFieldThetaData=thetaDataPerField{fi};
        else
             addThetaDataPerField
        end
        

        currFieldGoodSpikeTimes=spikeTimesInExpPerField{fi};
        
        currFieldThetaDataAllSpikeTimes=currFieldThetaData.allLapInFieldSpikeTimesInExp;
        allLapInFieldSpikeTraversalAvgSpeeds=currFieldThetaData.allLapInFieldSpikeTraversalAvgSpeeds;
        allLapInFieldSpikeTraversalSpeedStdDevs=currFieldThetaData.allLapInFieldSpikeTraversalSpeedStdDevs;
        
        currFieldThetaDataAllSpikeLaps=currFieldThetaData.allLapInFieldSpikeLapNums;
        
        currFieldThetaDataAllSpikeLaps=currFieldThetaDataAllSpikeLaps-min(currFieldThetaDataAllSpikeLaps)+1; %renorm at 1
        speedRangePerLap=currFieldThetaData.allLapInFieldSpeedRangePerTraversal;
        
        %goodLapNumIdxes=(speedRangePerLap<=maxSpeedStdDev);

        goodSpikeIdxesThetaData=zeros(size(currFieldThetaDataAllSpikeTimes));
        %lowSpeedVarAllSpikeIdxes=zeros(size(currFieldThetaDataAllSpikeTimes));
        
        for f=1:length(currFieldThetaDataAllSpikeTimes)
            currCandidateSpikeTime=currFieldThetaDataAllSpikeTimes(f);
            
            currCandidateSpikeLap=currFieldThetaDataAllSpikeLaps(f);
            
            %if curr spike is contained in list of good spike times
            if(min(abs(currFieldGoodSpikeTimes-currCandidateSpikeTime))<containsSpikeTol)
                goodSpikeIdxesThetaData(f)=1;
            end
            
            %{
            if(goodLapNumIdxes(currCandidateSpikeLap))
                lowSpeedVarAllSpikeIdxes(f)=1;
            end
            %}
        end
        
        currFieldThetaDataGoodSpikeIdxes=logical(goodSpikeIdxesThetaData);
        
        goodSpikeTraversalAvgSpeeds=allLapInFieldSpikeTraversalAvgSpeeds(currFieldThetaDataGoodSpikeIdxes);
        goodSpikeTraversalSpeedStdDevs=allLapInFieldSpikeTraversalSpeedStdDevs(currFieldThetaDataGoodSpikeIdxes);
        
        %lowSpeedVarSpikeIdxes=lowSpeedVarAllSpikeIdxes(currFieldThetaDataGoodSpikeIdxes);
        
        lowSpeedVarSpikeIdxes=goodSpikeTraversalSpeedStdDevs<=maxSpeedStdDev;
        
        if(excludeHighVarTraversals)
            goodSpikeTraversalAvgSpeeds(~lowSpeedVarSpikeIdxes)=NaN;
        end
        
        %goodSpikeTraversalAvgSpeeds(
        %goodSpikeTraversalSpeedRanges=(lowSpeedVarAllSpikeIdxes);
        
        posInfo=load(data.positionInfoForSessionSavePath);
        
        currFieldSpikePhases=spikePhasesPerField{fi};
        
        currFieldSpeedAlongTrackPerTimeMSec=posInfo.speedAlongTrackPerTimeMSec;
        positionTimeAxisSec=posInfo.positionTimeAxisSec;
        
        %currFieldSpikeSpeeds=interp1(positionTimeAxisSec,currFieldSpeedAlongTrackPerTimeMSec,currFieldSpikeTimesInExp);
        
        spikeAvgSpeedInTraversalsPerField{fi}=abs(goodSpikeTraversalAvgSpeeds(:));
        
        if(recalcThetaData)
            thetaDataPerField{fi}=currFieldThetaData;
            fieldNumWithinUnitPerField{fi}=fieldNumWithinUnit;
        end
        

         toc
    end

    %save('allUnitInFieldGoodSpeedSpikePhasesAndTimes.mat','spikeSpeedsPerField','thetaDataPerField','-append')
    if(~recalcThetaData)
         save('allUnitInFieldGoodSpeedSpikePhasesAndTimes.mat','spikeAvgSpeedInTraversalsPerField', '-append')
    else
        save('allUnitInFieldGoodSpeedSpikePhasesAndTimes.mat','spikeAvgSpeedInTraversalsPerField','thetaDataPerField', 'fieldNumWithinUnitPerField','-append')
    end

%else %recalc conditional
%    spikeDataPerField=load('allUnitInFieldGoodSpeedSpikePhasesAndTimes.mat');
end
%%
spikeDataPerField=load('allUnitInFieldGoodSpeedSpikePhasesAndTimes.mat');

%collect early phase vs avg speed for each field

totalFieldCount=spikeDataPerField.totalFieldCount;


earlyFieldTimeThresh=0.3;
earlyFieldTimeThresh=0.15;
earlyFieldTimeThresh=0.1;

earlyFieldTimeThresh=0.2;
%earlyFieldTimeThresh=0.15;

lowSpeedBounds=[0.2 0.4];

highSpeedBounds=[0.5 0.7];

lowSpeedBounds=[0.1 0.3];


lowSpeedBounds=[0 0.35];
highSpeedBounds=[0.35 0.7];
%highSpeedBounds=[0.6 0.8];


lowSpeedBounds=[0 0.3];
highSpeedBounds=[0.3 0.6];

lowSpeedBounds=[0 0.3];
highSpeedBounds=[0.3 0.5];

lowSpeedBounds=[0 0.3];
highSpeedBounds=[0.3 0.6];

zScoreSpeed=1;

if(zScoreSpeed)
    %lowSpeedBounds=[-2 -0.5];
    %highSpeedBounds=[0.5 2];
    
    lowSpeedBounds=[-3 -0.5];
    highSpeedBounds=[0.5 3];
end

%{
lowSpeedBounds=[0 0.2];
highSpeedBounds=[0.4 0.5];
%}


%{
lowSpeedBounds=[0.3 0.45];
highSpeedBounds=[0.45 0.6];
%}

%{
lowSpeedBounds=[0.2 0.3];
highSpeedBounds=[0.4 0.6];
%}

minNumSpikesInField=20;
fH=figure
circShift=-180;

fieldCountLowAndHighSpeed=0;
allFieldLowVsHighSpeedOffset=[];

phaseEdges=linspace(0,360,100);

usePeakPhase=1;
usePeakPhase=0;

minSpikesInEarlySpeedBin=3;

fEx=figure

setTightSubplots_Medium

startFi=1;
%startFi=492;

speedVsPhaseInfoPerField=cell(totalFieldCount,1);

showPlots=0;

for fi=startFi:totalFieldCount
    fi
    
    currFieldSpikeSpeeds=spikeDataPerField.spikeSpeedsPerField{fi};
    currFieldSpikeAvgSpeedInTraversals=spikeDataPerField.spikeAvgSpeedInTraversalsPerField{fi};
    
    highMinusLowPhaseDiffInfoCurrField={};
    
    unitData=load(spikeDataPerField.unitInfoPathPerField{fi});
    
    
    currFieldDirNum=spikeDataPerField.dirPerField{fi};
        
        if(currFieldDirNum==1)
            currFieldDirStr='Rightward';
        elseif(currFieldDirNum==2)
            currFieldDirStr='Leftward';
        end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %INSTANTANEOUS SPEED
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %currFieldSpikeAvgSpeedInTraversals=currFieldSpikeSpeeds;
    
    if(isempty(currFieldSpikeAvgSpeedInTraversals))
        continue
    end
    currFieldSpikeTimeFracInField=spikeDataPerField.spikeTimeFracInFieldPerField{fi};
    currFieldSpikePhases=spikeDataPerField.spikePhasesPerField{fi};
    
    
    
    currFieldEarlySpikeIdxes=currFieldSpikeTimeFracInField<earlyFieldTimeThresh;
    
    
    %currFieldSpikePhases=mod(currFieldSpikePhases+circShift,360);
    
    currFieldEarlySpikePhases=currFieldSpikePhases(currFieldEarlySpikeIdxes);
   
    currFieldEarlySpikeAvgTravSpeeds=currFieldSpikeAvgSpeedInTraversals(currFieldEarlySpikeIdxes);
    
    avgTravSpeedAcrossLaps=nanmean(unique(currFieldEarlySpikeAvgTravSpeeds));
    stdTravSpeedAcrossLaps=nanstd(unique(currFieldEarlySpikeAvgTravSpeeds));
    
    %zscore speed per field
    if(zScoreSpeed)
        currFieldSpikeAvgSpeedInTraversals=(currFieldSpikeAvgSpeedInTraversals-avgTravSpeedAcrossLaps)/stdTravSpeedAcrossLaps;
        currFieldEarlySpikeAvgTravSpeeds=(currFieldEarlySpikeAvgTravSpeeds-avgTravSpeedAcrossLaps)/stdTravSpeedAcrossLaps;
        
        lowSpeedBoundsRealSpeed=(lowSpeedBounds*stdTravSpeedAcrossLaps)+avgTravSpeedAcrossLaps;
        highSpeedBoundsRealSpeed=(highSpeedBounds*stdTravSpeedAcrossLaps)+avgTravSpeedAcrossLaps;
    else
        lowSpeedBoundsRealSpeed=lowSpeedBounds;
        highSpeedBoundsRealSpeed=highSpeedBounds;
    end
    currFieldLowSpeedSpikeIdxes=currFieldSpikeAvgSpeedInTraversals>=lowSpeedBounds(1) & currFieldSpikeAvgSpeedInTraversals<=lowSpeedBounds(2); 
    currFieldHighSpeedSpikeIdxes=currFieldSpikeAvgSpeedInTraversals>=highSpeedBounds(1) & currFieldSpikeAvgSpeedInTraversals<=highSpeedBounds(2); 

    currFieldLowSpeedSpikePhases=currFieldSpikePhases( currFieldLowSpeedSpikeIdxes);
    currFieldHighSpeedSpikePhases=currFieldSpikePhases( currFieldHighSpeedSpikeIdxes);
    
    currFieldLowSpeedSpikeTimeFracs=currFieldSpikeTimeFracInField( currFieldLowSpeedSpikeIdxes);
    currFieldHighSpeedSpikeTimeFracs=currFieldSpikeTimeFracInField( currFieldHighSpeedSpikeIdxes);
    
    
    uniqueTimeFracs=unique(currFieldSpikeTimeFracInField);
    uniqueTimeFracs(isnan(uniqueTimeFracs))=[];
    
    timeFracTol=0.001;
    minNumPhasesInTimeFrac=2;
    
    highMinusLowPhaseDiffPerTimeFrac=NaN(size(uniqueTimeFracs));
    for ti=1:length(uniqueTimeFracs)
        currTimeFrac=uniqueTimeFracs(ti);
        
        currTimeFracLowSpeedIdxes=(abs(currFieldLowSpeedSpikeTimeFracs-currTimeFrac)<timeFracTol);
        currTimeFracHighSpeedIdxes=(abs(currFieldHighSpeedSpikeTimeFracs-currTimeFrac)<timeFracTol);
        
        
        currTimeFracLowSpeedPhases=currFieldLowSpeedSpikePhases(currTimeFracLowSpeedIdxes);
        currTimeFracHighSpeedPhases=currFieldHighSpeedSpikePhases(currTimeFracHighSpeedIdxes);
        
      
        
        if(length(currTimeFracLowSpeedPhases)>=minNumPhasesInTimeFrac && length(currTimeFracHighSpeedPhases)>=minNumPhasesInTimeFrac)
              circMeanLowSpeedPhase=circMeanDeg(currTimeFracLowSpeedPhases);
                circMeanHighSpeedPhase=circMeanDeg(currTimeFracHighSpeedPhases);
                
            highMinusLowPhaseDiffPerTimeFrac(ti)=angdiffDeg([ circMeanLowSpeedPhase circMeanHighSpeedPhase]);
        end
        
    end

    thisFieldNumWithinUnit=spikeDataPerField.fieldNumWithinUnitPerField{fi};
    
    thisUnitStartEndM=unitData.(sprintf('%sFieldStartEndM',lower(currFieldDirStr)));
    
    thisFieldStartM=thisUnitStartEndM(2*(thisFieldNumWithinUnit-1)+1);
    thisFieldEndM=thisUnitStartEndM(2*(thisFieldNumWithinUnit-1)+2);
    
        
    allTraversalMeanSpeeds=spikeDataPerField.thetaDataPerField{fi}.allLapInFieldMeanSpeedPerTraversal;
    allTraversalSpeedStdDevs=spikeDataPerField.thetaDataPerField{fi}.allLapInFieldSpeedRangePerTraversal;
    
     lowTraversalSpeeds=allTraversalMeanSpeeds(allTraversalMeanSpeeds>=lowSpeedBoundsRealSpeed(1) & allTraversalMeanSpeeds<=lowSpeedBoundsRealSpeed(2));
      
     highTraversalSpeeds=allTraversalMeanSpeeds(allTraversalMeanSpeeds>=highSpeedBoundsRealSpeed(1) & allTraversalMeanSpeeds<=highSpeedBoundsRealSpeed(2));

    %highVsLowAbsSpeedDiff=nanmean(highSpeedBoundsRealSpeed)-nanmean(lowSpeedBoundsRealSpeed);
        highVsLowAbsSpeedDiff=nanmean(highTraversalSpeeds)-nanmean(lowTraversalSpeeds);
        
        meanHighSpeed=nanmean(highTraversalSpeeds);
        meanLowSpeed=nanmean(lowTraversalSpeeds);
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %collect offset info
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    currFieldLowSpeedEarlySpikePhases=currFieldSpikePhases(currFieldEarlySpikeIdxes & currFieldLowSpeedSpikeIdxes);
    currFieldHighSpeedEarlySpikePhases=currFieldSpikePhases(currFieldEarlySpikeIdxes & currFieldHighSpeedSpikeIdxes);
    
    
    currFieldAllEarlySpikePhases=currFieldSpikePhases(currFieldEarlySpikeIdxes);
    
    if(length(currFieldAllEarlySpikePhases)>=minSpikesInEarlySpeedBin)
        currFieldAllEarlyFieldMeanOffset=circMeanDeg(currFieldAllEarlySpikePhases);
        [~,currFieldAllEarlyFieldPeakPhase] = getCircKernelDistr(currFieldAllEarlySpikePhases,phaseEdges);
    else
        currFieldAllEarlyFieldMeanOffset=NaN;
        currFieldAllEarlyFieldPeakPhase=NaN;
    end

    %if(~(isempty(currFieldLowSpeedEarlySpikePhases) || isempty(currFieldHighSpeedEarlySpikePhases)))
        
    currFieldLowSpeedEarlyFieldMeanOffset=NaN;
        currFieldHighSpeedEarlyFieldMeanOffset=NaN;
        lowSpeedEarlyFieldPeakPhase=NaN;
        highSpeedEarlyFieldPeakPhase=NaN;
        
    if(length(currFieldLowSpeedEarlySpikePhases)>=minSpikesInEarlySpeedBin && length(currFieldHighSpeedEarlySpikePhases) >=minSpikesInEarlySpeedBin )
        currFieldLowSpeedEarlyFieldMeanOffset=circMeanDeg(currFieldLowSpeedEarlySpikePhases);
        currFieldHighSpeedEarlyFieldMeanOffset=circMeanDeg(currFieldHighSpeedEarlySpikePhases);

        [~,lowSpeedEarlyFieldPeakPhase] = getCircKernelDistr(currFieldLowSpeedEarlySpikePhases,phaseEdges);
        [~,highSpeedEarlyFieldPeakPhase] = getCircKernelDistr(currFieldHighSpeedEarlySpikePhases,phaseEdges);

        
        
        fieldCountLowAndHighSpeed=fieldCountLowAndHighSpeed+1;

        %{
        if(usePeakPhase)
            allFieldLowVsHighSpeedOffset=[allFieldLowVsHighSpeedOffset; [lowSpeedPeakPhase highSpeedPeakPhase]];
        else
           allFieldLowVsHighSpeedOffset=[allFieldLowVsHighSpeedOffset; [currFieldLowSpeedMeanOffset currFieldHighSpeedMeanOffset]];

        end
        %}

        if(showPlots)
            figure(fH)
            plot(currFieldEarlySpikeAvgTravSpeeds,currFieldEarlySpikePhases,'k.')
            hold on
             ylim([0 360])
            xlim([0 0.7])

            if(zScoreSpeed)
                xlim([-3 3])
            end
        end
    end
    
    speedEdges=linspace(0,0.7,30+1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %save curr field high speed vs low speed info
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     highMinusLowPhaseDiffInfoCurrField.currFieldLowSpeedSpikeTimeFracs=currFieldLowSpeedSpikeTimeFracs;
     highMinusLowPhaseDiffInfoCurrField.currFieldLowSpeedSpikePhases=currFieldLowSpeedSpikePhases;

      highMinusLowPhaseDiffInfoCurrField.currFieldHighSpeedSpikeTimeFracs=currFieldHighSpeedSpikeTimeFracs;
      highMinusLowPhaseDiffInfoCurrField.currFieldHighSpeedSpikePhases=currFieldHighSpeedSpikePhases;
      
      highMinusLowPhaseDiffInfoCurrField.highMinusLowPhaseDiffPerTimeFrac=highMinusLowPhaseDiffPerTimeFrac;
      highMinusLowPhaseDiffInfoCurrField.uniqueTimeFracs=uniqueTimeFracs;
      
       highMinusLowPhaseDiffInfoCurrField.lowSpeedBoundsZ=lowSpeedBounds;
      highMinusLowPhaseDiffInfoCurrField.highSpeedBoundsZ=highSpeedBounds;
      highMinusLowPhaseDiffInfoCurrField.lowSpeedBoundsRealSpeed=lowSpeedBoundsRealSpeed;
      highMinusLowPhaseDiffInfoCurrField.highSpeedBoundsRealSpeed=highSpeedBoundsRealSpeed;
      highMinusLowPhaseDiffInfoCurrField.meanLowSpeed=meanLowSpeed;
      highMinusLowPhaseDiffInfoCurrField.meanHighSpeed=meanHighSpeed;
      
      highMinusLowPhaseDiffInfoCurrField.allTraversalMeanSpeeds=allTraversalMeanSpeeds;
       highMinusLowPhaseDiffInfoCurrField.allTraversalSpeedStdDevs=allTraversalSpeedStdDevs;
      
      characteristicSpeed=nanmean(allTraversalMeanSpeeds);
      
      travSpeedStdDevAcrossLaps=nanstd(allTraversalMeanSpeeds);
      
      highMinusLowPhaseDiffInfoCurrField.characteristicSpeedMperSec=characteristicSpeed;
      highMinusLowPhaseDiffInfoCurrField.travSpeedStdDevAcrossLapsMperSec=travSpeedStdDevAcrossLaps;
      
      highMinusLowPhaseDiffInfoCurrField.currFieldLowSpeedEarlyFieldMeanOffset=currFieldLowSpeedEarlyFieldMeanOffset;
        highMinusLowPhaseDiffInfoCurrField.currFieldHighSpeedEarlyFieldMeanOffset=currFieldHighSpeedEarlyFieldMeanOffset;
        highMinusLowPhaseDiffInfoCurrField.lowSpeedEarlyFieldPeakPhase=lowSpeedEarlyFieldPeakPhase;
        highMinusLowPhaseDiffInfoCurrField.highSpeedEarlyFieldPeakPhase=highSpeedEarlyFieldPeakPhase;
        highMinusLowPhaseDiffInfoCurrField.minSpikesInEarlySpeedBin=minSpikesInEarlySpeedBin;
      
        highMinusLowPhaseDiffInfoCurrField.allEarlySpikePhases=currFieldAllEarlySpikePhases;
        highMinusLowPhaseDiffInfoCurrField.allEarlyFieldMeanOffset=currFieldAllEarlyFieldMeanOffset;
        highMinusLowPhaseDiffInfoCurrField.allEarlyFieldPeakPhase=currFieldAllEarlyFieldPeakPhase;
        highMinusLowPhaseDiffInfoCurrField.earlyFieldTimeThresh=earlyFieldTimeThresh;
  
     
    
      speedVsPhaseInfoPerField{fi}=highMinusLowPhaseDiffInfoCurrField;
      highMinusLowPhaseDiffInfoCurrField={};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %plots
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(showPlots)
        fEx=figure
        %subplot(10,10,fi)
        subplot(2,2,1)
        plot(currFieldLowSpeedSpikeTimeFracs,currFieldLowSpeedSpikePhases,'b.','MarkerSize',10)
        hold on
        plot(currFieldHighSpeedSpikeTimeFracs,currFieldHighSpeedSpikePhases,'r.','MarkerSize',10)


        xlabel('Time in field (frac)')
        xlim([0 1.2])
        ylim([0 360])
        ylabel('Theta phase')

        title({sprintf('Phase precession, low vs high speed for field %d',fi),sprintf('field bounds: %.2f m to %.2f m',thisFieldStartM,thisFieldEndM)})


        %plot(NaN,NaN,'r-')
        %plot(NaN,NaN,'b-')

        %lgd.FontSize=5;
        daspect([1 360 1])


        subplot(2,2,2)
   
        histogram(allTraversalMeanSpeeds,speedEdges)
        xlim([0 1])
        hold on

        %halfHeight=max(ylim)/2;
        plotHeight=max(ylim);

        vline(lowSpeedBoundsRealSpeed,'b--',3)

         vline(highSpeedBoundsRealSpeed,'r--',3)

            pR=plot(meanHighSpeed,plotHeight*0.8,'r.','MarkerSize',20)
            pB=plot(meanLowSpeed,plotHeight*0.8,'b.','MarkerSize',20)

        %ylim([0 50])
        xlabel('Average speed in traversal (m/s)')
        ylabel('Lap count')
        title({'Traversal speed distribution across laps',sprintf('High vs low absolute speed diff: %.2f m/s - %.2f m/s = %.2f m/s',meanHighSpeed,meanLowSpeed,highVsLowAbsSpeedDiff)})

          legend([pB pR], sprintf('traversal speed: %.2f to %.2f z',lowSpeedBounds(1),lowSpeedBounds(2)),sprintf('traversal speed: %.2f to %.2f z',highSpeedBounds(1),highSpeedBounds(2)),'Location','westoutside')
        legend boxoff

        subplot(2,2,[3])

        plot(uniqueTimeFracs, highMinusLowPhaseDiffPerTimeFrac,'r.','MarkerSize',30)

        hold on

        plot([0 1.2],[0 0],'k--')

        ylim([-180 180])
        xlim([0 1.2])

            xlabel('Time in field (frac)')
        ylabel('(avg high speed phase) - (avg low speed phase) (deg)')
        title({'High - low speed phase difference across field'})
        daspect([1 360 1])

        subplot(2,2,[4])

        try
            imshow(unitData.unitInfo.(sprintf('saveImgPath%s',currFieldDirStr)))
        end

        setFigFontTo(16)
        maxFig

        saveas(gcf,fullfile(exampleFieldsDir,sprintf('lowVsHighSpeedPhasePrecession_Field%d.png',fi)))
        close(fEx)

        %if(justExamplePlots)
        %    continue
        %end
        drawnow

    end



        %{
        plot(1,currFieldLowSpeedMeanOffset,'k.')
        hold on
        plot(2,currFieldHighSpeedMeanOffset,'k.')
        %}

        %analyze continuous relationship between speed and early phase



        %{
        if(length(currFieldEarlySpikePhases)>=minNumSpikesInField)
            %[rho,p,slopeDegPerXunit, offsetDeg]=getCircCorrCoeff(currFieldEarlySpikeAvgTravSpeeds,currFieldEarlySpikePhases,fH)
            [m,b,R,p]=getLinearFit(currFieldEarlySpikeAvgTravSpeeds,currFieldEarlySpikePhases)

        end

        %}

          
    
end
%%

save('allUnitInFieldGoodSpeedSpikePhasesAndTimes.mat','speedVsPhaseInfoPerField', '-append')
%close all
phaseBins=edgesToBins(phaseEdges);

lowToHighDiffs=NaN(fieldCountLowAndHighSpeed,1);
for fi=1:fieldCountLowAndHighSpeed
    currFieldLowToHighDiff=angdiffDeg([allFieldLowVsHighSpeedOffset(fi,1),allFieldLowVsHighSpeedOffset(fi,2)]);
   %currFieldLowToHighDiff=diff([allFieldLowVsHighSpeedOffset(fi,1),allFieldLowVsHighSpeedOffset(fi,2)]);

    lowToHighDiffs(fi)=currFieldLowToHighDiff;
end

figure;
minDiff=-120;
maxDiff=120;
[Nc,diffEdges]=histcounts(lowToHighDiffs,linspace(minDiff,maxDiff,50))
h=histogram(lowToHighDiffs,linspace(minDiff,maxDiff,50))
%bar(edgesToBins(diffEdges),Nc)
hold on
try
ylim([0 max(Nc)])
xlim([minDiff maxDiff])
end

zH=vline(0,'k--',3)

mH=vline(nanmean(lowToHighDiffs),'r--',5)

legend([zH mH],{'0 difference','mean difference'})

legend boxoff

xlabel('high - low speed phase offset per field (degrees)')
ylabel('Field count')

[pVal,h,stats] = signrank(lowToHighDiffs);

zVal=stats.zval;


title(sprintf('Wilcoxon signed rank test, z=%.2f, p=%d, n=%d fields',zVal,pVal,fieldCountLowAndHighSpeed))

[lowSpeedPhaseOffsetDistr,lowSpeedAllFieldPeakPhase]=getCircKernelDistr(allFieldLowVsHighSpeedOffset(:,1),phaseEdges);
[highSpeedPhaseOffsetDistr,highSpeedAllFieldPeakPhase]=getCircKernelDistr(allFieldLowVsHighSpeedOffset(:,2),phaseEdges);
setFigFontTo(18)

saveas(gcf,'thetaOffsetVsSpeedPerIndividualFieldDiffDistr.png')

figure; 
pLow=plot(phaseBins,lowSpeedPhaseOffsetDistr,'b-','LineWidth',3);
hold on;
pHigh=plot(phaseBins,highSpeedPhaseOffsetDistr,'r-','LineWidth',3);

legend('low speed (-2 < z < -0.5)', 'high speed (0.5 < z < 2)','Location','northwest')
legend boxoff
xlabel('Circ-mean theta phase in early field (deg)')
ylabel('Probability')
axis tight
setFigFontTo(18)
title(sprintf('Theta offset vs speed per individual field, n=%d fields',fieldCountLowAndHighSpeed))

saveas(gcf,'thetaOffsetVsSpeedPerIndividualFieldDistr.png')


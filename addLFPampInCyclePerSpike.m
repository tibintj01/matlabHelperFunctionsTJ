%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load unit struct file names
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; clc
processedDataDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData';
unitDataDir=fullfile(processedDataDir,'unitSpikeInfoStructs');

saveLFPcycleInfoDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData/lfpCycleInfos';
saveRawLFPInfoDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData/rawLFPstructs';

filePaths=getFilePathsRegex(unitDataDir,'*mat');

spikeDataPerField=load('allUnitInFieldGoodSpeedSpikePhasesAndTimes.mat')

recalcLFPperSpike=0;
recalcLFPperSpike=1;


showPlots=0;

justGammaAppend=1;
if(showPlots)
    justGammaAppend=0;
end

if(recalcLFPperSpike)

    spikeTimesInExpPerField=spikeDataPerField.spikeTimesInExpPerField;
    spikePhasesPerField =spikeDataPerField.spikePhasesPerField;
    unitInfoPathPerField=spikeDataPerField.unitInfoPathPerField;

    totalFieldCount=spikeDataPerField.totalFieldCount;

    startFieldIdx=1;
     fH=figure;

     allFieldPhaseVsNormInCycleLFPVsZscoreLFPvsGammaDurCycleFrac=[];
     spikeLFPvaluesNormInCyclePerField=cell(size(spikeTimesInExpPerField));
        spikeLFPzScoreValuesPerField=cell(size(spikeTimesInExpPerField));
        spikeGammaDurationsThetaFracPerField=cell(size(spikeTimesInExpPerField));
        spikeGammaDurationsSecPerField=cell(size(spikeTimesInExpPerField));
        
    for fi=startFieldIdx:totalFieldCount
        %for fi=76:length(filePaths)
        tic
           fi
        currFilePath=unitInfoPathPerField{fi};
        currFileName=getFileNameFromPath(currFilePath);
        fileBaseName=currFileName(1:(end-4));

        dataStruct=load(currFilePath);
        data=dataStruct.unitInfo;

        currFieldSpikeTimesInExp=spikeTimesInExpPerField{fi};
        currFieldSpikePhases=spikePhasesPerField{fi};

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %load pyramidal layer lfp
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        maxRippleChThisShank=data.maxRippleChThisShank;
        currSessionName=data.sessionName;
       rawLFPfilePath=getRegexFilePath(saveRawLFPInfoDir,sprintf('*_%s_fullLFPstruct.mat', currSessionName));


       currSessPyrLayerRawLFPs=load(rawLFPfilePath);
       
  

       currUnitShankNum=data.unitShankNum;

       lfpShankNums=currSessPyrLayerRawLFPs.pyrLayerChPerShank(:,1);
       lfpShankRefChNums=currSessPyrLayerRawLFPs.pyrLayerChPerShank(:,2);

       currRefCh=lfpShankRefChNums(currUnitShankNum==lfpShankNums);

       timeAxisRawLFP=currSessPyrLayerRawLFPs.timeAxis;
       currRefChRawLFP=currSessPyrLayerRawLFPs.rawLFPpyrLayerChs.(sprintf('ch%d',currRefCh));
       
       %lowFreqGamma=30;
            lowFreqGamma=20;
    %highFreqGamma=120;
        highFreqGamma=150;
        minGammAmpZ=-1;
        maxGammAmpZ=1;
        
    cycleInfoFileName=fullfile(saveLFPcycleInfoDir,sprintf('%s_Ch%d_%d-%d_minmax.mat',currSessionName,currRefCh,lowFreqGamma,highFreqGamma));
    currUnitGammaCycleInfo=load(cycleInfoFileName);

    gammaMinTimes=currUnitGammaCycleInfo.cycleMinTimes;
    gammaMaxTimes=currUnitGammaCycleInfo.cycleMaxTimes;
    
    gammaAmps=currUnitGammaCycleInfo.cycleMaxAmp-currUnitGammaCycleInfo.cycleMinAmp;
    gammaAmpZscore=zscoreLFP(gammaAmps);
    
    goodGammaCycleIdxes=gammaAmpZscore>minGammAmpZ & gammaAmpZscore<maxGammAmpZ;
    
    gammaMaxTimes(~goodGammaCycleIdxes)=NaN;
    gammaMinTimes(~goodGammaCycleIdxes)=NaN;
    
       
       %currChGammaCycleInfo=
       
       currRefChRawZscoreLFP=currSessPyrLayerRawLFPs.rawZscoreLFPpyrLayerChs.(sprintf('ch%d',currRefCh));

       lowFreq=6;
       highFreq=12;
        cycleInfoFileName=fullfile(saveLFPcycleInfoDir,sprintf('%s_Ch%d_%d-%d_minmax.mat',currSessionName,maxRippleChThisShank,lowFreq,highFreq));
        currUnitCycleInfo=load(cycleInfoFileName);

        minTimes=currUnitCycleInfo.cycleMinTimes;
        maxTimes=currUnitCycleInfo.cycleMaxTimes;


        numCycles=length(currUnitCycleInfo.cycleMaxTimes);

        %[spikeAssignMax spikeAssignMin spikeAssignZero] = assignSpikesToCycles2017(currFieldSpikeTimesInExp, minTimes, maxTimes);

        currFieldSpikeLFPvaluesNormInCycle=NaN(size(currFieldSpikeTimesInExp));
        currFieldSpikeLFPzScoreValues=NaN(size(currFieldSpikeTimesInExp));
        currFieldSpikeGammaDurationsSec=NaN(size(currFieldSpikeTimesInExp));
        currFieldSpikeGammaDurationsThetaFrac=NaN(size(currFieldSpikeTimesInExp));
        
        dtLFP=median(diff(timeAxisRawLFP));
        smoothWindSec=0.005; %5 msec
        
        smoothWindIdx=round(smoothWindSec/dtLFP);
        for ci=1:(numCycles-1)
            currCycleStartTime=maxTimes(ci);
            currCycleEndTime=maxTimes(ci+1);
            
            currCycleDurationSec=currCycleEndTime-currCycleStartTime;

            currCycleSpikeTimeIdxes=currFieldSpikeTimesInExp<=currCycleEndTime & currFieldSpikeTimesInExp>=currCycleStartTime;
            spikeTimesInCurrCycle=currFieldSpikeTimesInExp(currCycleSpikeTimeIdxes);

            gammaMaxTimesCurrCycle=gammaMaxTimes(gammaMaxTimes<=currCycleEndTime & gammaMaxTimes>=currCycleStartTime);
            gammaMaxIdxesAfter=find(gammaMaxTimes>currCycleEndTime);
            gammaMaxIdxesBefore=find(gammaMaxTimes<currCycleStartTime);
            
            if(~isempty(gammaMaxIdxesAfter))
                firstGammaMaxTimeAfter=gammaMaxTimes(gammaMaxIdxesAfter(1));
            else
                firstGammaMaxTimeAfter=NaN;
            end
            if(~isempty(gammaMaxIdxesBefore))
             lastGammaMaxTimeBefore=gammaMaxTimes(gammaMaxIdxesBefore(end));
            else
                lastGammaMaxTimeBefore=NaN;
            end
            
           gammaMaxTimesCurrCyclePlusMinusOne=[lastGammaMaxTimeBefore  ;gammaMaxTimesCurrCycle(:); firstGammaMaxTimeAfter];
           
            %gammaMaxTimesCurrCycle=gammaMinTimes(gammaMinTimes<=currCycleEndTime & gammaMinTimes>=currCycleStartTime);

            if(isempty(spikeTimesInCurrCycle))
                continue
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %get LFP value (norm in cycle) per spike
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if(~justGammaAppend)
                currCycleLFPtimes=timeAxisRawLFP(timeAxisRawLFP>=currCycleStartTime & timeAxisRawLFP<=currCycleEndTime);
                currCycleLFPvalues=currRefChRawLFP(timeAxisRawLFP>=currCycleStartTime & timeAxisRawLFP<=currCycleEndTime);
            end
            
            if(showPlots)
                plot(1000*(currCycleLFPtimes-min(currCycleLFPtimes)),smooth(currCycleLFPvalues,smoothWindIdx))

                 hold on
                 plot(1000*(spikeTimesInCurrCycle-min(currCycleLFPtimes)),zeros(size(spikeTimesInCurrCycle)),'k*')

                 vline(1000*(gammaMaxTimesCurrCyclePlusMinusOne-min(currCycleLFPtimes)),'r')
            end
             
            
             numGammaCycles=length(gammaMaxTimesCurrCyclePlusMinusOne);
             spikeGammaCycleIdxes=floor(interp1(gammaMaxTimesCurrCyclePlusMinusOne,1:numGammaCycles,spikeTimesInCurrCycle));
             
             gammaDurations=diff(gammaMaxTimesCurrCyclePlusMinusOne);
             
             if(min(gammaDurations)<0)
                 disp('')
                 
             end
             
             currCycleSpikeGammaDurationsSec=gammaDurations(spikeGammaCycleIdxes);
             
             currCycleSpikeGammaDurationsThetaFrac=currCycleSpikeGammaDurationsSec/currCycleDurationSec;
             
            
             
             if(~justGammaAppend)
                currCycleRawLFPvalues=interp1(timeAxisRawLFP,currRefChRawLFP,currCycleLFPtimes);



                currCycleAmp=range(currCycleRawLFPvalues);
                currCycleMin=min(currCycleRawLFPvalues);
                currCycleMax=max(currCycleRawLFPvalues);



                currCycleSpikeTimesRawLFPvalues=interp1(timeAxisRawLFP,currRefChRawLFP,spikeTimesInCurrCycle);

                currCycleSpikeTimesRawLFPvaluesNormInCycle=(currCycleSpikeTimesRawLFPvalues-currCycleMin)/currCycleAmp;

                currCycleSpikeTimesRawLFPzScoreValues=interp1(timeAxisRawLFP,currRefChRawZscoreLFP,spikeTimesInCurrCycle);


               currFieldSpikeLFPzScoreValues(currCycleSpikeTimeIdxes)=currCycleSpikeTimesRawLFPzScoreValues;
                currFieldSpikeLFPvaluesNormInCycle(currCycleSpikeTimeIdxes)=currCycleSpikeTimesRawLFPvaluesNormInCycle;
            
             end
            currFieldSpikeGammaDurationsThetaFrac(currCycleSpikeTimeIdxes)=currCycleSpikeGammaDurationsThetaFrac;
            currFieldSpikeGammaDurationsSec(currCycleSpikeTimeIdxes)=currCycleSpikeGammaDurationsSec;
            
            if(showPlots)
                close all
            end
            
        end

        if(~justGammaAppend)
            spikeLFPvaluesNormInCyclePerField{fi}=currFieldSpikeLFPvaluesNormInCycle(:);
            spikeLFPzScoreValuesPerField{fi}=currFieldSpikeLFPzScoreValues(:);
            allFieldPhaseVsNormInCycleLFPVsZscoreLFPvsGammaDurCycleFrac=[allFieldPhaseVsNormInCycleLFPVsZscoreLFPvsGammaDurCycleFrac; [currFieldSpikePhases(:), currFieldSpikeLFPvaluesNormInCycle(:), currFieldSpikeLFPzScoreValues(:), currFieldSpikeGammaDurationsThetaFrac(:)]];

        end
        spikeGammaDurationsThetaFracPerField{fi}=currFieldSpikeGammaDurationsThetaFrac(:);
        spikeGammaDurationsSecPerField{fi}=currFieldSpikeGammaDurationsSec(:);
        

         toc
    end

    if(~justGammaAppend)
        save('allUnitInFieldGoodSpeedSpikePhasesAndTimes.mat','allFieldPhaseVsNormInCycleLFPVsZscoreLFPvsGammaDurCycleFrac','spikeLFPvaluesNormInCyclePerField','spikeLFPzScoreValuesPerField','spikeGammaDurationsThetaFracPerField','spikeGammaDurationsSecPerField','-append')
    else
      save('allUnitInFieldGoodSpeedSpikePhasesAndTimes.mat','spikeGammaDurationsThetaFracPerField','spikeGammaDurationsSecPerField','-append')

    end
else %recalc conditional
    load('allUnitInFieldGoodSpeedSpikePhasesAndTimes.mat')
end
%%
%lfpValueEdges=linspace(-3,3,201);

actualTime=0;
actualTime=1;
if(actualTime)
    lfpValueEdges=linspace(0.004,0.04,31);
    %lfpValueEdges=linspace(0,120,31);
else
    lfpValueEdges=linspace(0.05,0.2,31);
    lfpValueEdges=linspace(0.04,0.4,71);
end

phaseEdges=linspace(0,360,31);

%phaseEdges=linspace(0,360,11);
fH=figure; 
circShiftDeg=0;

lfpCycleNorm1Zscore2Gamma3=3;


%allFieldPhaseVsLFPvaluesCircShift=mod(allFieldPhaseVsNormInCycleLFPVsZscoreLFPvsGammaDurCycleFrac(:,1)+circShiftDeg,360);
allFieldPhaseVsLFPvaluesCircShift=vertcat(spikePhasesPerField{:});

if(actualTime)
    allFieldGammaDurationsThetaFrac=vertcat(spikeGammaDurationsSecPerField{:});

else
    allFieldGammaDurationsThetaFrac=vertcat(spikeGammaDurationsThetaFracPerField{:});
end

getJointDistrGivenX(allFieldPhaseVsLFPvaluesCircShift,allFieldGammaDurationsThetaFrac,phaseEdges, lfpValueEdges,fH,0)
xlim([0 360])
%ylim([-2.25 2.25])

%caxis([0 2e-4])
%caxis([0 3.5e-4])
xlabel('Theta phase of in-field spike (deg)')

if(lfpCycleNorm1Zscore2Gamma3==1)
    ylabel('LFP value (norm. within cycle)')
elseif(lfpCycleNorm1Zscore2Gamma3==2)
    ylabel('LFP value (z-score across session)')
else
    ylabel('Gamma cycle duration (frac of theta cycle)')
end

title(sprintf('Local gamma cycle duration vs theta phase of spikes, n=%d fields',totalFieldCount))
setFigFontTo(18)

%plot(allFieldPhaseVsLFPvalues(:,1),allFieldPhaseVsLFPvalues(:,2),'k.')
     hold on;
     drawnow
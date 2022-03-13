%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load unit struct file names
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; clc
processedDataDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData';
unitDataDir=fullfile(processedDataDir,'unitSpikeInfoStructs');

saveLFPcycleInfoDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData/lfpCycleInfos';
saveRawLFPInfoDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData/rawLFPstructs';

filePaths=getFilePathsRegex(unitDataDir,'*mat');

spikeDataPerField=load('allUnitInFieldGoodSpeedSpikePhasesAndTimes.mat');

recalcLFPperSpike=0;
recalcLFPperSpike=1;

minDataPerBin=2;

minDataPerBin=1;

cycleBufferTime= 50/1000; %sec

periSpikeBufferTimeSec=30/1000;
%periSpikeBufferTimeSec=15/1000;

periSpikeBufferTimeSec=100;
periSpikeBufferTimeSec=10/1000;

periSpikeBufferTimeSec=5/1000;

showPlots=0;
showPlotsPerField=0;

%justGammaAppend=1;
justGammaAppend=0;
if(showPlots || showPlotsPerField)
    justGammaAppend=0;
end

if(recalcLFPperSpike)

    spikeTimesInExpPerField=spikeDataPerField.spikeTimesInExpPerField;
    spikePhasesPerField =spikeDataPerField.spikePhasesPerField;
    unitInfoPathPerField=spikeDataPerField.unitInfoPathPerField;

    totalFieldCount=spikeDataPerField.totalFieldCount;
    
    timeInFieldPerCyclePerField=spikeDataPerField.timeInFieldPerCyclePerField;
    

    startFieldIdx=1;
     fH=figure;

     allFieldPhaseVsNormInCycleLFPVsZscoreLFPvsGammaDurCycleFrac=[];
     spikeLFPvaluesNormInCyclePerField=cell(size(spikeTimesInExpPerField));
        spikeLFPzScoreValuesPerField=cell(size(spikeTimesInExpPerField));
        spikeGammaDurationsThetaFracPerField=cell(size(spikeTimesInExpPerField));
        spikeGammaDurationsSecPerField=cell(size(spikeTimesInExpPerField));
        
        periSpikeGammaWTacrossThetaPerField=cell(totalFieldCount,1);
        
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
        
        currFieldTimeInFieldFracPerCycle=timeInFieldPerCyclePerField{fi};

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
       lfpFS=median(diff(timeAxisRawLFP));
       currRefChRawLFP=currSessPyrLayerRawLFPs.rawLFPpyrLayerChs.(sprintf('ch%d',currRefCh));
       
       %{
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
    
       %}
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
        
        wtCommonGammaFreqAxis=(25:120)';
            %wtCommonThetaPhaseAxis=0:5:355;
             wtCommonThetaPhaseAxis=1:1:360;
            
            numFreqBins=length(wtCommonGammaFreqAxis);
            
        periSpikeWTavgAcrossThetaCycles=zeros(numFreqBins,length(wtCommonThetaPhaseAxis));
        numWTs=0;
        dataCountGrid=zeros(numFreqBins,length(wtCommonThetaPhaseAxis));
     
        
        nonSpikeWTpowerPerFreqPerCycle=NaN(numFreqBins,numCycles);
        
        for ci=1:(numCycles-1)
            currCycleStartTime=maxTimes(ci);
            currCycleEndTime=maxTimes(ci+1);
            
            currCycleDurationSec=currCycleEndTime-currCycleStartTime;

            currCycleSpikeTimeIdxes=currFieldSpikeTimesInExp<=currCycleEndTime & currFieldSpikeTimesInExp>=currCycleStartTime;
            spikeTimesInCurrCycle=currFieldSpikeTimesInExp(currCycleSpikeTimeIdxes);
            
        

            %currFieldSection=
            %all theta cycles part of section
            if(isempty(spikeTimesInCurrCycle))
                continue
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %get LFP value (norm in cycle) per spike
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if(~justGammaAppend)
                currCycleLFPtimes=timeAxisRawLFP(timeAxisRawLFP>=currCycleStartTime & timeAxisRawLFP<=currCycleEndTime);
                currCycleLFPvalues=currRefChRawLFP(timeAxisRawLFP>=currCycleStartTime & timeAxisRawLFP<=currCycleEndTime);
                
                currCycleWithBufferLFPtimes=timeAxisRawLFP(timeAxisRawLFP>=currCycleStartTime-cycleBufferTime & timeAxisRawLFP<=currCycleEndTime+cycleBufferTime);
                currCycleWithBufferLFPvalues=currRefChRawLFP(timeAxisRawLFP>=currCycleStartTime-cycleBufferTime & timeAxisRawLFP<=currCycleEndTime+cycleBufferTime);
            end
            
            originalCycleIdxes=(currCycleWithBufferLFPtimes>=currCycleStartTime & currCycleWithBufferLFPtimes<=currCycleEndTime);
            
            if(showPlots)
                plot(1000*(currCycleLFPtimes-min(currCycleLFPtimes)),smooth(currCycleLFPvalues,smoothWindIdx))

                 hold on
                 plot(1000*(spikeTimesInCurrCycle-min(currCycleLFPtimes)),zeros(size(spikeTimesInCurrCycle)),'k*')

                 vline(1000*(gammaMaxTimesCurrCyclePlusMinusOne-min(currCycleLFPtimes)),'r')
            end
             
            
         
             
            
             
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
           
            
            
           
          
 
            
            
            %[thisCycleWT,freqs] = cwt(currCycleLFPvalues,1./lfpFS);
            [thisCycleWT,freqs] = cwt(currCycleWithBufferLFPvalues,1./lfpFS);
            
            thisCycleWT=thisCycleWT(:,originalCycleIdxes);
            
            
                
            dispFreqIdxes=freqs<=100;
            
            
            currCycleLFPphases=360*(currCycleLFPtimes-currCycleStartTime)/currCycleDurationSec;
            
            thisCycleWTinterp=interp2(currCycleLFPphases,freqs,thisCycleWT,wtCommonThetaPhaseAxis,wtCommonGammaFreqAxis);
            
            periSpikeIdxes=zeros(size(wtCommonThetaPhaseAxis));
            
             for si=1:length(spikeTimesInCurrCycle)
                currSpikeTime=spikeTimesInCurrCycle(si);
                
                periSpikeMinTime=max(currSpikeTime-periSpikeBufferTimeSec,currCycleStartTime);
                periSpikeMaxTime=min(currCycleEndTime, currSpikeTime+periSpikeBufferTimeSec);

                periSpikeMinPhase=360*(periSpikeMinTime-currCycleStartTime)/currCycleDurationSec;
                periSpikeMaxPhase=360*(periSpikeMaxTime-currCycleStartTime)/currCycleDurationSec;

                periSpikeIdxes=periSpikeIdxes |  (wtCommonThetaPhaseAxis<=periSpikeMaxPhase & wtCommonThetaPhaseAxis>=periSpikeMinPhase); %accumulate regions around spikes
            
             end
            
            nonSpikeIdxes=~periSpikeIdxes;
            
           totalPower=sum(abs(thisCycleWTinterp(~isnan(thisCycleWTinterp(:,periSpikeIdxes)))));
           
               nonSpikeWTpowerPerFreqPerCycle(:,ci)=nanmean(abs(thisCycleWTinterp(:,nonSpikeIdxes)),2)/totalPower;
            
            
            
            thisCycleWTinterp(:,~periSpikeIdxes)=NaN;
            
             numWTs=numWTs+1;
             
             dataCountGrid=dataCountGrid+(~isnan(thisCycleWTinterp));
             thisCycleWTinterp(isnan(thisCycleWTinterp))=0;
             
             
             
             periSpikeWTavgAcrossThetaCycles=periSpikeWTavgAcrossThetaCycles+(abs(thisCycleWTinterp)/totalPower);
             
             
           
             %periSpikeWTavgAcrossThetaCycles=periSpikeWTavgAcrossThetaCycles+real(thisCycleWTinterp);

             
             

             
           
             
            highCountIdxes=dataCountGrid>10;
             
        end
        
        meanNonSpikePowerPerFreq=nanmean(nonSpikeWTpowerPerFreqPerCycle,2);
        stdNonSpikePowerPerFreq=nanstd(nonSpikeWTpowerPerFreqPerCycle,[],2);
        
        %periSpikeWTavgAcrossThetaCycles=(periSpikeWTavgAcrossThetaCycles-meanNonSpikePowerPerFreq)./stdNonSpikePowerPerFreq;
        
        periSpikeWTavgAcrossThetaCycles=periSpikeWTavgAcrossThetaCycles/nansum(periSpikeWTavgAcrossThetaCycles(:));
        %periSpikeWTavgAcrossThetaCycles=(periSpikeWTavgAcrossThetaCycles-nanmean(periSpikeWTavgAcrossThetaCycles(:)))/nanstd(periSpikeWTavgAcrossThetaCycles(:));
        periSpikeWTavgAcrossThetaCycles(dataCountGrid<minDataPerBin)=NaN;
         if(showPlotsPerField)
                
                %figure;
            



            %omarPcolor(currCycleLFPtimes-currCycleStartTime,freqs(dispFreqIdxes),abs(thisCycleWT(dispFreqIdxes,:)))
            omarPcolor(wtCommonThetaPhaseAxis,wtCommonGammaFreqAxis,periSpikeWTavgAcrossThetaCycles./dataCountGrid)

            maxC=prctile(periSpikeWTavgAcrossThetaCycles(highCountIdxes)./dataCountGrid(highCountIdxes),99.99);
            colormap(jet)
            cb=colorbar
            ylabel(cb,'wavelet power')

            if(~isnan(maxC))
             %caxis([0 maxC])
            end

            autoArrangeFigures()
            close all
               
         end
        
         periSpikeGammaWTacrossThetaPerField{fi}=periSpikeWTavgAcrossThetaCycles;
         


        
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
    
    periSpikeGammaWTacrossThetaPerField_JustAroundSpikes=periSpikeGammaWTacrossThetaPerField;
    save('allUnitInFieldGoodSpeedSpikePhasesAndTimes.mat','periSpikeGammaWTacrossThetaPerField_JustAroundSpikes','wtCommonThetaPhaseAxis','wtCommonGammaFreqAxis','periSpikeBufferTimeSec', '-append')
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

%getJointDistrGivenX(allFieldPhaseVsLFPvaluesCircShift,allFieldGammaDurationsThetaFrac,phaseEdges, lfpValueEdges,fH,0)
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
     
     gammaWaveletAvgAnalysis
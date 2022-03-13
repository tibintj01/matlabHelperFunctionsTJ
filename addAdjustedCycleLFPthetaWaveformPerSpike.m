%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load unit struct file names
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;  clc

%NOT RELOADING LARGE FILE, TAKE OUT TO RELOAD
%clear all;
clearvars -except spikeDataPerField


processedDataDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData';
unitDataDir=fullfile(processedDataDir,'unitSpikeInfoStructs');

saveLFPcycleInfoDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData/lfpCycleInfos';
saveRawLFPInfoDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData/rawLFPstructs';

useOmarAsymAdjustment=1;

if(useOmarAsymAdjustment)
    adjustedCycleLFPInfoDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData/lfpOmarAsymAdjustedCycleInfos';

else
       adjustedCycleLFPInfoDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData/lfpAdjustedCycleInfos';

end

filePaths=getFilePathsRegex(unitDataDir,'*mat');

tic
disp('loading spike data....')
if(~exist('spikeDataPerField','var'))
    spikeDataPerField=load('allUnitInFieldGoodSpeedSpikePhasesAndTimes.mat');
end
toc
recalcLFPperSpike=0;
recalcLFPperSpike=1;

showPlots=0;
%showPlots=1;


if(recalcLFPperSpike)

    spikeTimesInExpPerField=spikeDataPerField.spikeTimesInExpPerField;
    spikePhasesPerField =spikeDataPerField.spikePhasesPerField;
    unitInfoPathPerField=spikeDataPerField.unitInfoPathPerField;

    totalFieldCount=spikeDataPerField.totalFieldCount;

    startFieldIdx=1;
     fH=figure;
     
    adjustedThetaWaveformLFPzPerSpikePerField=cell(size(spikeTimesInExpPerField));
    
    adjustedThetaSpikePhasesPerField=cell(size(spikeTimesInExpPerField));
        
    for fi=startFieldIdx:totalFieldCount
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
        %load pyramidal layer lfp cycle info
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        maxRippleChThisShank=data.maxRippleChThisShank;
        currSessionName=data.sessionName;
       rawLFPfilePath=getRegexFilePath(saveRawLFPInfoDir,sprintf('*_%s_fullLFPstruct.mat', currSessionName));
        adjustedLFPfilePath=getRegexFilePath(adjustedCycleLFPInfoDir,sprintf('*_%s_adjCycleLFPstruct.mat', currSessionName));

       currSessPyrLayerRawLFPs=load(rawLFPfilePath);
       currSessPyrLayerAdjustedCycleInfos=load(adjustedLFPfilePath);

       currUnitShankNum=data.unitShankNum;

       lfpShankNums=currSessPyrLayerRawLFPs.pyrLayerChPerShank(:,1);
       lfpShankRefChNums=currSessPyrLayerRawLFPs.pyrLayerChPerShank(:,2);

       currRefCh=lfpShankRefChNums(currUnitShankNum==lfpShankNums);

       timeAxisRawLFP=currSessPyrLayerRawLFPs.timeAxis;
       currRefChRawLFP=currSessPyrLayerRawLFPs.rawLFPpyrLayerChs.(sprintf('ch%d',currRefCh));
       currRefChRawZscoreLFP=currSessPyrLayerRawLFPs.rawZscoreLFPpyrLayerChs.(sprintf('ch%d',currRefCh));
       currRefChAdjustedCycleInfo=currSessPyrLayerAdjustedCycleInfos.adjustedCycleInfoPerPyrCh.(sprintf('ch%d',currRefCh));

       %{
       lowFreq=6;
       highFreq=12;
        cycleInfoFileName=fullfile(saveLFPcycleInfoDir,sprintf('%s_Ch%d_%d-%d_minmax.mat',currSessionName,maxRippleChThisShank,lowFreq,highFreq));
        %}
        %currUnitCycleInfo=load(cycleInfoFileName);

        if(useOmarAsymAdjustment)
            minTimes=currRefChAdjustedCycleInfo.fastMidTimes;

            adjustedMaxTimes=currRefChAdjustedCycleInfo.fastStartTimes;

            numCycles=length(currRefChAdjustedCycleInfo.fastStartTimes);
        else
            minTimes=currRefChAdjustedCycleInfo.cycleMinTimes;

            adjustedMaxTimes=currRefChAdjustedCycleInfo.closestFastMaxTimes;

            numCycles=length(currRefChAdjustedCycleInfo.closestFastMaxTimes);
        end

        currFieldThetaWaveformLFPzPerSpike=cell(size(currFieldSpikeTimesInExp));
        
        currFieldAdjustedThetaSpikePhases=NaN(size(currFieldSpikeTimesInExp));

        
        dtLFP=median(diff(timeAxisRawLFP));
        smoothWindSec=0.005; %5 msec
        
        smoothWindIdx=round(smoothWindSec/dtLFP);
        
        for ci=1:(numCycles-1)
            currCycleStartTime=adjustedMaxTimes(ci);
            currCycleEndTime=adjustedMaxTimes(ci+1);
            
            currCycleDurationSec=currCycleEndTime-currCycleStartTime;

            currCycleSpikeTimeIdxBool=currFieldSpikeTimesInExp<=currCycleEndTime & currFieldSpikeTimesInExp>=currCycleStartTime;
            spikeTimesInCurrCycle=currFieldSpikeTimesInExp(currCycleSpikeTimeIdxBool);

            currCycleSpikeTimeIdxes=find(currCycleSpikeTimeIdxBool);

            if(isempty(spikeTimesInCurrCycle))
                continue
            end
            
            adjustedSpikePhasesInCurrCycle=360*(spikeTimesInCurrCycle-currCycleStartTime)/currCycleDurationSec;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %get LFP waveform per spike
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            

                currCycleLFPtimes=timeAxisRawLFP(timeAxisRawLFP>=currCycleStartTime & timeAxisRawLFP<=currCycleEndTime);
                %currCycleLFPvalues=currRefChRawLFP(timeAxisRawLFP>=currCycleStartTime & timeAxisRawLFP<=currCycleEndTime);
                currCycleLFPvalues=currRefChRawZscoreLFP(timeAxisRawLFP>=currCycleStartTime & timeAxisRawLFP<=currCycleEndTime);


            
                currCycleLFPtimeAxis=currCycleLFPtimes-min(currCycleLFPtimes);
                
                for si=1:length(currCycleSpikeTimeIdxes)
                    currSpikeIdx=currCycleSpikeTimeIdxes(si);
                    currFieldThetaWaveformLFPzPerSpike{currSpikeIdx}=[currCycleLFPtimeAxis(:), currCycleLFPvalues(:)];
                    
                    currFieldAdjustedThetaSpikePhases(currSpikeIdx)=adjustedSpikePhasesInCurrCycle(si);
                end
                
            if(showPlots)
                plot(1000*currFieldThetaWaveformLFPzPerSpike{currSpikeIdx}(:,1),smooth(currFieldThetaWaveformLFPzPerSpike{currSpikeIdx}(:,2),smoothWindIdx))
                drawnow
                close all
            end
            
        end

        adjustedThetaWaveformLFPzPerSpikePerField{fi}=currFieldThetaWaveformLFPzPerSpike;    
        adjustedThetaSpikePhasesPerField{fi}=currFieldAdjustedThetaSpikePhases;

         toc
    end
    
    if(useOmarAsymAdjustment)
        omarAsymAdjustedThetaWaveformLFPzPerSpikePerField=adjustedThetaWaveformLFPzPerSpikePerField;
        omarAsymAdjustedThetaSpikePhasesPerField=adjustedThetaSpikePhasesPerField;
        save('allUnitInFieldGoodSpeedSpikePhasesAndTimes.mat','omarAsymAdjustedThetaWaveformLFPzPerSpikePerField','omarAsymAdjustedThetaSpikePhasesPerField','-append')

    else
      save('allUnitInFieldGoodSpeedSpikePhasesAndTimes.mat','adjustedThetaWaveformLFPzPerSpikePerField','adjustedThetaSpikePhasesPerField','-append')
    end
else %recalc conditional
    load('allUnitInFieldGoodSpeedSpikePhasesAndTimes.mat')
end
%%
%lfp
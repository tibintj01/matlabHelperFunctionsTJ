%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load unit struct file names
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; clc
processedDataDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData';
unitDataDir=fullfile(processedDataDir,'unitSpikeInfoStructs');

saveLFPcycleInfoDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData/lfpCycleInfos';
saveRawLFPInfoDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData/rawLFPstructs';
saveRawLFPInfoDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData/lfpAdjustedCycleInfos';


filePaths=getFilePathsRegex(unitDataDir,'*mat');

spikeDataPerField=load('allUnitInFieldGoodSpeedSpikePhasesAndTimes.mat');

recalcLFPperSpike=0;
recalcLFPperSpike=1;

showPlots=0;


if(recalcLFPperSpike)

    spikeTimesInExpPerField=spikeDataPerField.spikeTimesInExpPerField;
    spikePhasesPerField =spikeDataPerField.spikePhasesPerField;
    unitInfoPathPerField=spikeDataPerField.unitInfoPathPerField;

    totalFieldCount=spikeDataPerField.totalFieldCount;

    startFieldIdx=1;
     fH=figure;
     
    thetaWaveformLFPzPerSpikePerField=cell(size(spikeTimesInExpPerField));
        
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
       
       
       currRefChRawZscoreLFP=currSessPyrLayerRawLFPs.rawZscoreLFPpyrLayerChs.(sprintf('ch%d',currRefCh));

       lowFreq=6;
       highFreq=12;
        cycleInfoFileName=fullfile(saveLFPcycleInfoDir,sprintf('%s_Ch%d_%d-%d_minmax.mat',currSessionName,maxRippleChThisShank,lowFreq,highFreq));
        currUnitCycleInfo=load(cycleInfoFileName);

        minTimes=currUnitCycleInfo.cycleMinTimes;
        maxTimes=currUnitCycleInfo.cycleMaxTimes;


        numCycles=length(currUnitCycleInfo.cycleMaxTimes);

        currFieldThetaWaveformLFPzPerSpike=cell(size(currFieldSpikeTimesInExp));

        
        dtLFP=median(diff(timeAxisRawLFP));
        smoothWindSec=0.005; %5 msec
        
        smoothWindIdx=round(smoothWindSec/dtLFP);
        
        for ci=1:(numCycles-1)
            currCycleStartTime=maxTimes(ci);
            currCycleEndTime=maxTimes(ci+1);
            
            currCycleDurationSec=currCycleEndTime-currCycleStartTime;

            currCycleSpikeTimeIdxBool=currFieldSpikeTimesInExp<=currCycleEndTime & currFieldSpikeTimesInExp>=currCycleStartTime;
            spikeTimesInCurrCycle=currFieldSpikeTimesInExp(currCycleSpikeTimeIdxBool);

            currCycleSpikeTimeIdxes=find(currCycleSpikeTimeIdxBool);

            if(isempty(spikeTimesInCurrCycle))
                continue
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %get LFP value (norm in cycle) per spike
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            

                currCycleLFPtimes=timeAxisRawLFP(timeAxisRawLFP>=currCycleStartTime & timeAxisRawLFP<=currCycleEndTime);
                %currCycleLFPvalues=currRefChRawLFP(timeAxisRawLFP>=currCycleStartTime & timeAxisRawLFP<=currCycleEndTime);
                currCycleLFPvalues=currRefChRawZscoreLFP(timeAxisRawLFP>=currCycleStartTime & timeAxisRawLFP<=currCycleEndTime);


            
                currCycleLFPtimeAxis=currCycleLFPtimes-min(currCycleLFPtimes);
                
                for si=1:length(currCycleSpikeTimeIdxes);
                    currSpikeIdx=currCycleSpikeTimeIdxes(si);
                    currFieldThetaWaveformLFPzPerSpike{currSpikeIdx}=[currCycleLFPtimeAxis(:), currCycleLFPvalues(:)];
                end
                
            if(showPlots)
                plot(1000*currFieldThetaWaveformLFPzPerSpike{currSpikeIdx}(:,1),smooth(currFieldThetaWaveformLFPzPerSpike{currSpikeIdx}(:,2),smoothWindIdx))
                drawnow
                close all
            end
            
        end

        thetaWaveformLFPzPerSpikePerField{fi}=currFieldThetaWaveformLFPzPerSpike;        

         toc
    end
      save('allUnitInFieldGoodSpeedSpikePhasesAndTimes.mat','thetaWaveformLFPzPerSpikePerField','-append')

else %recalc conditional
    load('allUnitInFieldGoodSpeedSpikePhasesAndTimes.mat')
end
%%
%lfp
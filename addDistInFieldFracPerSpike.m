%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load unit struct file names
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;  clc

%NOT RELOADING LARGE FILE, TAKE OUT TO RELOAD
clear all;
%clearvars -except spikeDataPerField

processedDataDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData';
unitDataDir=fullfile(processedDataDir,'unitSpikeInfoStructs');

filePaths=getFilePathsRegex(unitDataDir,'*mat');

tic
disp('loading spike data....')
if(~exist('spikeDataPerField','var'))
    spikeDataPerField=load('allUnitInFieldGoodSpeedSpikePhasesAndTimes.mat');
end
toc
recalcLFPperSpike=0;
%recalcLFPperSpike=1;

%showPlots=0;
showPlots=1;

if(recalcLFPperSpike)

    spikeTimesInExpPerField=spikeDataPerField.spikeTimesInExpPerField;
    spikePhasesPerField =spikeDataPerField.spikePhasesPerField;
    unitInfoPathPerField=spikeDataPerField.unitInfoPathPerField;
   fieldNumWithinUnitPerField=spikeDataPerField.fieldNumWithinUnitPerField;
   dirPerField=spikeDataPerField.dirPerField;


    totalFieldCount=spikeDataPerField.totalFieldCount;

    startFieldIdx=1;
     fH=figure;
     
    distInFieldPerSpikePerField=cell(size(spikeTimesInExpPerField));
        
    for fi=startFieldIdx:totalFieldCount
         tic
         
           fi
        currFilePath=unitInfoPathPerField{fi};
        currFileName=getFileNameFromPath(currFilePath);
        fileBaseName=currFileName(1:(end-4));

        dataStruct=load(currFilePath);


        posInfo=load(dataStruct.unitInfo.positionInfoForSessionSavePath);
        posPerTimeM=posInfo.posAlongTrackPerTimeM;
        positionTimeAxisSec=posInfo.positionTimeAxisSec;
        
        currFieldNumWithinUnit=fieldNumWithinUnitPerField{fi};
        currFieldDir=dirPerField{fi};
        
        currFieldStartIdx=currFieldNumWithinUnit*2-1;
        if(currFieldDir==1) %rightward
            currFieldBoundsM=dataStruct.rightwardFieldStartEndM(currFieldStartIdx:currFieldStartIdx+1);
        elseif(currFieldDir==2)
            currFieldBoundsM=dataStruct.leftwardFieldStartEndM(currFieldStartIdx:currFieldStartIdx+1);
        end
        
        currFieldWidthM=abs(diff(currFieldBoundsM));
        
        currFieldSpikeTimesInExp=spikeTimesInExpPerField{fi};
        
        currFieldPosPerSpike=interp1(positionTimeAxisSec,posPerTimeM,currFieldSpikeTimesInExp);
        currFieldSpikePhases=spikePhasesPerField{fi};
        
        if(currFieldDir==1) %rightward
              currFieldDistInFieldPerSpike=(currFieldPosPerSpike-currFieldBoundsM(1))/currFieldWidthM;
        elseif(currFieldDir==2)
            currFieldDistInFieldPerSpike=-(currFieldPosPerSpike-currFieldBoundsM(1))/currFieldWidthM;
        end
        
        distInFieldPerSpikePerField{fi}=currFieldDistInFieldPerSpike;

         toc
    end

    save('allUnitInFieldGoodSpeedSpikePhasesAndTimes.mat','distInFieldPerSpikePerField','-append')

else %recalc conditional
    load('allUnitInFieldGoodSpeedSpikePhasesAndTimes.mat')
end




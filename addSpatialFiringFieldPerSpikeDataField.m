%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load unit struct file names
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clc;
clearvars -except spikeDataPerField

processedDataDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData';
unitDataDir=fullfile(processedDataDir,'unitSpikeInfoStructs');

%spikeDataPerField=load('allUnitInFieldGoodSpeedSpikePhasesAndTimes.mat');

totalFieldCount=spikeDataPerField.totalFieldCount;
firingRatePerPosPerField=cell(totalFieldCount,1);

for fi=1:totalFieldCount
    fi
    currFieldUnitData=load(spikeDataPerField.unitInfoPathPerField{fi});
    
    currFieldDirection=spikeDataPerField.dirPerField{fi}; %1= rightward, 2 = leftward
    
    posInfo=load(currFieldUnitData.unitInfo.positionInfoForSessionSavePath);
    
    fieldSpikeTimes=spikeDataPerField.spikeTimesInExpPerField{fi};
    
    if(currFieldDirection==1) %rightward
       % currFieldFiringRatePerPos=currFieldUnitData.unitInfo.firingRatePerPositionRight;
        
        [positionBins,firingRatePerPositionRight, ~,~] ...
               = getSpikeStatsPerPosition(fieldSpikeTimes,posInfo.positionTimeAxisSec,posInfo.posAlongTrackPerTimeM,posInfo.speedAlongTrackPerTimeMSec);
    
           currFieldFiringRatePerPos=firingRatePerPositionRight;
           currFieldPosBins=positionBins;
    elseif(currFieldDirection==2) %leftward
        %currFieldFiringRatePerPos=currFieldUnitData.unitInfo.firingRatePerPositionLeft;
        [positionBins,~,firingRatePerPositionLeft,~] ...
               = getSpikeStatsPerPosition(fieldSpikeTimes,posInfo.positionTimeAxisSec,posInfo.posAlongTrackPerTimeM,posInfo.speedAlongTrackPerTimeMSec);
    
           currFieldFiringRatePerPos=firingRatePerPositionLeft;
           currFieldPosBins=positionBins;
           
    end
    %currFieldPosBins=currFieldUnitData.unitInfo.positionBins;
    
    currFieldFiringRatePerPos=[currFieldPosBins(:), currFieldFiringRatePerPos(:)];
    
    
    
    firingRatePerPosPerField{fi}=currFieldFiringRatePerPos;
end

save('allUnitInFieldGoodSpeedSpikePhasesAndTimes.mat','firingRatePerPosPerField','-append')

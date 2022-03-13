%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load unit struct file names
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;  clc

%NOT RELOADING LARGE FILE, TAKE OUT TO RELOAD
%clear all;
clearvars -except spikeDataPerField

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
recalcLFPperSpike=1;

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
     
    accPerSpikePerField=cell(size(spikeTimesInExpPerField));
        
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
        
        currFieldDir=dirPerField{fi};
        
        
        currFieldSpikeTimesInExp=spikeTimesInExpPerField{fi};
 
                derivModelOrder=2; %local quadratic fit
         
                dt=median(diff(positionTimeAxisSec));
                       derivativeEstWind=round(2/dt); %seconds to points
                        %derivativeEstWind=round(1/dt); %seconds to points

                
                if(currFieldDir==1)
                     speedPerTimeMperSec=movingslope(posPerTimeM,derivativeEstWind,derivModelOrder,dt);
                elseif(currFieldDir==2)
                    speedPerTimeMperSec=-movingslope(posPerTimeM,derivativeEstWind,derivModelOrder,dt); %leftward speed
                end
               
              
               accPerTimeMperSecSquared=movingslope(speedPerTimeMperSec,derivativeEstWind,derivModelOrder,dt);
       
               
        currFieldAccPerSpike=interp1(positionTimeAxisSec,accPerTimeMperSecSquared,currFieldSpikeTimesInExp);
  
      
        
        accPerSpikePerField{fi}=currFieldAccPerSpike;

         toc
    end

    save('allUnitInFieldGoodSpeedSpikePhasesAndTimes.mat','accPerSpikePerField','-append')

else %recalc conditional
    load('allUnitInFieldGoodSpeedSpikePhasesAndTimes.mat')
end




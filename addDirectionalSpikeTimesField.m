%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load unit struct file names
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; clc
processedDataDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData';
unitDataDir=fullfile(processedDataDir,'unitSpikeInfoStructs');


filePaths=getFilePathsRegex(unitDataDir,'*mat');

startFile=1;
for fi=startFile:length(filePaths)
       fi
    currFilePath=filePaths{fi};
    currFileName=getFileNameFromPath(currFilePath);
    fileBaseName=currFileName(1:(end-4));
    
    dataStruct=load(currFilePath);
    data=dataStruct.unitInfo;
     
    
    for di=1:2

         currDirSpikeTimes=data.spikeTimes(data.spikeTimeDirectionAssignment==di);

         if(di==2)
            leftSpikeTimes=currDirSpikeTimes;
            save(currFilePath,'leftSpikeTimes','-append')
         elseif(di==1)
            rightSpikeTimes=currDirSpikeTimes;
            save(currFilePath,'rightSpikeTimes','-append')
         end
        
        
    end

end
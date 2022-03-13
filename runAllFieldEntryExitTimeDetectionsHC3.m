close all; clear all; clc

processedDataDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData';
markedBoundsSaveDir='./hc3Images/phasePrecessionBounds';

touchDir(markedBoundsSaveDir)

unitDataDir=fullfile(processedDataDir,'unitSpikeInfoStructs');

filePaths=getFilePathsRegex(unitDataDir,'*mat');

for fi=1:length(filePaths)

    currFilePath=filePaths{fi};
    fi/length(filePaths)
    
    %if(strContainsCircSessionName(currFilePath))
       % try
            getFieldEntryExitTimesPerLapHC3(currFilePath);
      %{
	  catch ME
             disp('SKIPPING%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
            disp(ME.message)
            disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        end
	%}
    %end
end

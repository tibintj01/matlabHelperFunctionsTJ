close all; clear all; clc
dataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles';

filePaths=getFilePathsRegex(dataDir,'*Info.mat');

for fi=1:length(filePaths)
%for fi=1:1
    %for fi=98:length(filePaths)
    currFilePath=filePaths{fi};
    fi/length(filePaths)
    
    %if(strContainsCircSessionName(currFilePath))
       % try
            getFieldEntryExitTimesPerLap(currFilePath);
      %{
	  catch ME
             disp('SKIPPING%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
            disp(ME.message)
            disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        end
	%}
    %end
end

function [refUnitDataLeftward,refUnitDataRightward] = getRefUnitData(sessionName)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
 dataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles/';
 saveFileName=sprintf('%sLeftAndRightRefUnits.mat',sessionName); 
 if(exist(saveFileName,'file'))
     savedData=load(saveFileName);
     refUnitDataLeftward=savedData.refUnitDataLeftward;
     refUnitDataRightward=savedData.refUnitDataRightward;
     
     return
 else
     
 fileNameRegex=sprintf('%s_6-12Theta*Info.mat',sessionName);
    
    sampleDataFiles=getRegexFilePaths(dataDir,fileNameRegex);
       refUnitDataLeftward=NaN;
      refUnitDataRightward=NaN;
    for i=1:length(sampleDataFiles)
        currUnitData=load(sampleDataFiles{i});
        
        if(isfield(currUnitData,'manualFieldStartsMrightward') && isfield(currUnitData.directionSpecificStats,'rightward'))
            refUnitDataRightward=currUnitData;
            break
        end
    end
    
     
    for i=1:length(sampleDataFiles)
        currUnitData=load(sampleDataFiles{i});
        
        if(isfield(currUnitData,'manualFieldStartsMleftward')&& isfield(currUnitData.directionSpecificStats,'leftward'))
            refUnitDataLeftward=currUnitData;
            break
        end
       
    end
    
    
   save(saveFileName,'refUnitDataLeftward','refUnitDataRightward')
    

    end


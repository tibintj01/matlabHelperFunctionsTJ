close all; clear all; clc

saveDir='./hc3ProcessedData/unitSpikeInfoStructs';


unitInfoPaths=getFilePathsRegex(saveDir,'*mat');

for fi=1:length(unitInfoPaths)
    currUnitData=load(unitInfoPaths{fi});
    
    currUnitData.unitInfo
    
    
end
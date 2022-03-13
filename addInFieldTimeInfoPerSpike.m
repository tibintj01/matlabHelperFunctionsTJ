function [] = addInFieldTimeInfoPerSpike(currFilePath)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    cycleInfoDataDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData/lfpCycleInfos';

    unitStruct=load(currFilePath);
    
    refChCycleInfo=load(unitStruct.refChCycleInfoPath);
    
    unitData=unitStruct.unitInfo;
    sessionName=unitData.sessionName;
    
    
    
    %refChCycleInfoPath=getFilePathsRegex(cycleInfoDataDir,sprintf('%s_Ch%d*',sessionName,refChNum));
    
    %save(currFilePath,'-append')
    
    %updatedStruct=load(currFilePath)
    disp('')
    
    
    
    

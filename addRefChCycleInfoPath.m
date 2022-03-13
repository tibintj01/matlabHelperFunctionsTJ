function [] = addRefChCycleInfoPath(currFilePath)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    cycleInfoDataDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData/lfpCycleInfos';

    unitStruct=load(currFilePath);
    unitData=unitStruct.unitInfo;
    
    
    refChNum=unitData.maxRippleChThisShank;
    sessionName=unitData.sessionName;
    
    refChCycleInfoPath=getFilePathsRegex(cycleInfoDataDir,sprintf('%s_Ch%d_*.mat',sessionName,refChNum));
    
    save(currFilePath,'refChCycleInfoPath','-append')
    
    %updatedStruct=load(currFilePath)
    %disp('')
    
    
    

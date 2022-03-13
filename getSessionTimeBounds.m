function [minTime, maxTime,numLapsPerDir]=getSessionTimeBounds(sessionName)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    %sessionNames={'Achilles_10252013','Achilles_11012013','Buddy_06272013', 'Cicero_09102014','Gatsby_08282013'};

    dataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles/'
    fileNameRegex=sprintf('%s_6-12Theta*Info.mat',sessionName);
    
    sampleDataFiles=getRegexFilePaths(dataDir,fileNameRegex);
    
    currData=load(sampleDataFiles{1});
    
    minTime=min(currData.positionTimeAxis);
    maxTime=max(currData.positionTimeAxis);
    numLapsPerDir=currData.numLapsPerDir;
    


function [mostLockingThetaCh] = getRefThetaChNumForSession(fileName)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    [currSessionNum,currSessName]=getSessionNum(fileName);

%sessionNames={'Achilles_10252013','Achilles_11012013','Buddy_06272013','Cicero_09012014','Cicero_09102014','Cicero_09172014','Gatsby_08282013'};


    refThetaChannels=[103 63 80 NaN 127 NaN 4];
    
    mostLockingThetaCh=refThetaChannels(currSessionNum);
    disp(sprintf('retrieved reference theta channel for %s, Ch %d', currSessName,mostLockingThetaCh))
    
    


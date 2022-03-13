
lfpMatDir='/nfs/turbo/lsa-ojahmed/tibin/processedHumanData'
%MG63
%sessionIdx=6;

%convertNS5sToMatForSessionNum(sessionIdx)
%createConcatenatedLFPforSessionNum(sessionIdx,lfpMatDir)

%MG29
%RIHE1
%sessionIdx=13;
%convertNS5sToMatForSessionNum(sessionIdx)
createConcatenatedLFPforSessionNum(sessionIdx,lfpMatDir)

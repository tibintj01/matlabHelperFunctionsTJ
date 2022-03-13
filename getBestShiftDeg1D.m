function[bestShiftDeg,shiftedPhases]=getBestShiftDeg1D(localSpikePhasePerTime,localNormTimeInFieldPerTime,earlyTimeThresh)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    earlyThresh=0.2;
    earlyThresh=0.05;
     earlyThresh=earlyTimeThresh;
    %earlyThresh=0.25;
    %lateThresh=0.8
    
    lateThresh=1-earlyThresh;
    %earlyThresh=0.25;
    %lateThresh=0.75;
    earlyIdxes= localNormTimeInFieldPerTime(:) < earlyThresh;
     lateIdxes=localNormTimeInFieldPerTime(:) > lateThresh;
         %earlyIdxes=localNormDistInFieldPerTime(:)<earlyThresh;
         %lateIdxes=localNormDistInFieldPerTime(:)>lateThresh;

    
    maxShift=360;
    bestShiftDeg=0;
    numShifts=maxShift;
    linearAvgEarlyPhase=NaN(numShifts,1);
    
    %for i=1:numShifts
    for i=1:maxShift
        shiftedPhases=mod(localSpikePhasePerTime-i,360);
        %currAvg=nanmean(shiftedPhases(earlyIdxes));
        %currAvg=nanmedian(shiftedPhases(earlyIdxes));
        currAvgEarly=nanmean(shiftedPhases(earlyIdxes));
        currAvgLate=nanmean(shiftedPhases(lateIdxes));
         %currAvgEarly=nanmedian(shiftedPhases(earlyIdxes));
        %currAvgLate=nanmedian(shiftedPhases(lateIdxes));
        %linearAvgEarlyPhase(i)=currAvg;
         linearAvgEarlyPhase(i)=currAvgEarly-currAvgLate;
    end
    %linearAvgEarlyPhase(linearAvgEarlyPhase>maxShift)=NaN;
    [~,bestShiftDeg]=max(linearAvgEarlyPhase);
    
    shiftedPhases=mod(localSpikePhasePerTime-bestShiftDeg,360);
    
end


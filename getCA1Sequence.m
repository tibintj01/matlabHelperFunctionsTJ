function [ca1ActivationSequence,peakActivityStrengths]=getCA1Sequence(ca3ActivitySeries,postPreConn)
%
    ca1ActivitySeries=postPreConn*ca3ActivitySeries'; %(nCA1 x nCA3) x (nCA3 x numSteps)
    [peakActivityStrengths,mostActiveCellPerStep]=max(ca1ActivitySeries,[],1);
    ca1ActivationSequence=mostActiveCellPerStep;
end


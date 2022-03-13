function[bestShiftDeg,shiftedPhases]=getBestShiftDegZeroPhase(localSpikePhasePerTime,localNormTimeInFieldPerTime)
%finds most probable zerophase and shifts until this is maximal

    zeroTimeIdxes=localNormTimeInFieldPerTime==0;
    
    zeroTimePhases=localSpikePhasePerTime(zeroTimeIdxes);
    
    
    maxShift=360;
    bestShiftDeg=0;
    numShifts=maxShift;
    mostProbPhase=NaN(numShifts,1);
    linearRPerShift=NaN(numShifts,1);
    phaseEdges=linspace(0,360,361);
    
    %for i=1:numShifts
    for i=1:maxShift
        %shiftedPhases=mod(zeroTimePhases-i,360);
        shiftedPhases=mod(localSpikePhasePerTime-i,360);
        %currAvg=nanmean(shiftedPhases(earlyIdxes));
        %currAvg=nanmedian(shiftedPhases(earlyIdxes));
        %{
        currPdist=getCircKernelDistr(shiftedPhases,phaseEdges);
        [~,maxProbID]=max(currPdist);
        mostProbPhase(i)=maxProbID;
        %}
      
        
         [m,b,R]=getLinearFit(localNormTimeInFieldPerTime,shiftedPhases);
         linearRPerShift(i)=R;
    end
    %linearAvgEarlyPhase(linearAvgEarlyPhase>maxShift)=NaN;
    %[~,bestShiftDeg]=max(mostProbPhase);
    [~,bestShiftDeg]=min(linearRPerShift);
    
    
    shiftedPhases=mod(localSpikePhasePerTime-bestShiftDeg,360);
    
    
    
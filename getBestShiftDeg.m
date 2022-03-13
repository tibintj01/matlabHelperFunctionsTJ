function[bestShiftDeg,bestShiftedPhases]=getBestShiftDeg(localSpikePhasePerTime,localNormDistInFieldPerTime,localNormTimeInFieldPerTime)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    earlyThresh=0.2;
        lateThresh=0.8;
        
           % earlyThresh=0.1;
        %lateThresh=0.9;
        
        %earlyThresh=0.3;
        %lateThresh=0.7;
       %earlyThresh=0.3;
    %earlyThresh=0.25;
    
    %includeLateRestr=0;
    includeLateRestr=1;
    
   

    %lateThresh=1;
    
    %earlyThresh=0.25;
    %lateThresh=0.75;
    earlyIdxes=localNormDistInFieldPerTime(:)<=earlyThresh & localNormTimeInFieldPerTime(:) <= earlyThresh;
         lateIdxes=localNormDistInFieldPerTime(:)>=lateThresh & localNormTimeInFieldPerTime(:) >= lateThresh;
         
         %earlyIdxes=localNormDistInFieldPerTime(:)<=earlyThresh | localNormTimeInFieldPerTime(:) <= earlyThresh;
         %lateIdxes=localNormDistInFieldPerTime(:)>=lateThresh | localNormTimeInFieldPerTime(:) >= lateThresh;
         %earlyIdxes=localNormDistInFieldPerTime(:)<earlyThresh;
         %lateIdxes=localNormDistInFieldPerTime(:)>lateThresh;

    
    maxShift=360;
    bestShiftDeg=0;
    numShifts=maxShift;
    linearAvgEarlyPhase=NaN(numShifts,1);
    
    %for i=1:numShifts
    for i=1:maxShift
        shiftedPhases=mod(localSpikePhasePerTime-i,360);
        
        
  
        currAvgEarly=nanmean(shiftedPhases(earlyIdxes));
        currAvgLate=nanmean(shiftedPhases(lateIdxes));
        
     
         if(includeLateRestr)
            linearAvgEarlyPhase(i)=currAvgEarly-currAvgLate;
         else
             linearAvgEarlyPhase(i)=currAvgEarly;
         end
    end
    %linearAvgEarlyPhase(linearAvgEarlyPhase>maxShift)=NaN;
    [~,bestShiftDeg]=max(linearAvgEarlyPhase);
    

    bestShiftedPhases=mod(localSpikePhasePerTime-bestShiftDeg,360);


function [bestPhaseShift]=getBestLinearCorrShift(phasesDeg,positions)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %circularly shift phases until error from model linear precession is
        %minimal
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        degShifts=1:360;
        rhos=NaN(size(degShifts));
        totalVars=NaN(size(degShifts));
        for si=1:length(degShifts)
            currShift=degShifts(si);
            shiftedPhases=mod(phasesDeg+currShift,360);
            nonNaNids=~isnan(shiftedPhases);
            totalVars(si)=var(shiftedPhases(nonNaNids));
            currRho=corr(positions(nonNaNids),shiftedPhases(nonNaNids));
            rhos(si)=currRho;
        end
        [~,maxRhoID]=max(abs(rhos));
        %[~,minVarID]=min(totalVars);
        bestPhaseShift=degShifts(maxRhoID);

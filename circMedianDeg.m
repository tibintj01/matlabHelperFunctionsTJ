function [angMedian]=circMedianDeg(spikePhasesDeg,dim)
	if(~exist('dim','var'))
        %angMedian=mod(rad2ang(circ_median(ang2rad(spikePhasesDeg(:)))),360);
        notNaNIdxes=~isnan(spikePhasesDeg);
        
        nonNanPhasesDeg=spikePhasesDeg(notNaNIdxes);
          angMedian=mod(rad2ang(circ_median(ang2rad(nonNanPhasesDeg(:)))),360);
    else
        angMedian=mod(rad2ang(circ_median(ang2rad(spikePhasesDeg),[],dim)),360);
    end

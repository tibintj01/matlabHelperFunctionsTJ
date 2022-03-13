function [angMean]=circMeanDeg(spikePhasesDeg,dim)
	if(~exist('dim','var'))
        angMean=mod(rad2ang(circ_mean(ang2rad(spikePhasesDeg(:)))),360);
    else
        angMean=mod(rad2ang(circ_mean(ang2rad(spikePhasesDeg),[],dim)),360);
    end

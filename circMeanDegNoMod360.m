function [angMean]=circMeanDegNoMod360(spikePhasesDeg,dim)
	if(~exist('dim','var'))
        angMean=rad2ang(circ_mean(ang2rad(spikePhasesDeg(:))));
    else
        angMean=rad2ang(circ_mean(ang2rad(spikePhasesDeg),[],dim));
    end

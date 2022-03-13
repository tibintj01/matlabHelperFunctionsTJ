function [p,z]=rayleighZDeg(spikePhasesDeg)

	[p,z]=circ_rtest(ang2rad(spikePhasesDeg(:)));

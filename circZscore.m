function [circZ]=circZscoreDeg(angs)
	mu=circMeanDeg(angs);
	r=circResLength(ang2rad(angs));

	circStd=(sqrt(-2*log(r)))

	circZ=(angs-mu)/circStd;

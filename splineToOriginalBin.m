function [lowResBin]=splineToOriginalBin(highResBin)
	splineDt=0.1;

	lowResBin = (highResBin - 1/splineDt + 1) .* splineDt;

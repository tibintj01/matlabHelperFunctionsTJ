function [highResBin]=splineToOriginalBin(lowResBin)
	splineDt=0.1;

	highResBin=round(lowResBin/splineDt-1+1/splineDt);
	%lowResBin = (highResBin - 1/splineDt + 1) .* splineDt;

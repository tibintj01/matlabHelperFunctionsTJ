function [shiftedVals]=alignFirstNvals(vals,trueVals,firstN)
	firstN
	numNegShifts=round(firstN/2);
	numPosShifts=round(firstN/2);

	absErrorPosShift=zeros(1+numPosShifts,1);
	absErrorNegShift=zeros(1+numNegShifts,1);

	%compare values for which derivative is negative
	valsDiff=[0 diff(vals)];
	for i=0:numPosShifts
		shiftedVals=shiftVals(vals,i);
		absErrorPosShift(i)=getAbsError(vals(valsDiff<0),shiftedVals(valsDiff<0));
	end
	for i=1:numNegShifts
		shiftedVals=shiftVals(vals,-i);
		absErrorNegShift(i)=getAbsError(vals(valsDiff<0),shiftedVals(valsDiff<0));
	end


	allErrors=[absErrorNegShift(:);absErrorPosShift(:)];

	[bestError,bestShift]=min(allErrors);

	shiftedVals=shiftVals(vals,bestShift);


	

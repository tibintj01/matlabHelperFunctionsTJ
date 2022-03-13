function [centerIdx]=getCenterIdx(values)
	numer=0;
	denom=sum(values);

	for i=1:length(values)
		numer=numer+i*values(i);
	end

	centerIdx=round(numer/denom);

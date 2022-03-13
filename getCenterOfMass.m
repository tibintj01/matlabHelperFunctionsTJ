function [cmX,cmVal]=getCenterOfMass(values)

	numer=0;
	for i=1:length(values)
		numer=numer+values(i)*i;
	end

	denom=0;

	for i=1:length(values)
		denom=denom+values(i);
	end

	cmX=(numer/denom);
	cmVal=mean(values);

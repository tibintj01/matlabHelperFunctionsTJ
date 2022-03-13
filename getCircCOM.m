function [COM]=getCircCOM(x)

	COM=circMean(x,(1:length(x))'/length(x));

	return
	%{
	numer=0;

	for i=1:length(x)
		numer=numer+x(i)*i;
	end

	denom=0;
	for i=1:length(x)
		denom=denom+x(i);
	end
	COM=numer/denom;
	%}

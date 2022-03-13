function [Abar]=subtractMeanCol(A)
	n=size(A,2);
        Amean=A*ones(n,1)*ones(n,1)'/n;
	Abar=A-Amean;
	

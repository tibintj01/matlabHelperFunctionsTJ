function [y] = getUVecComb(timeAxis,coeffs)
	%combines first k left singular vectors using coeffs and returns
	%coeffs must be k-length vector

	freqBandName='Alpha';

	svdData=load(sprintf('singleCycle%sSVD.mat',freqBandName));
    	U=svdData.U;

	k=length(coeffs);
	Uk=U(:,1:k);
	y=Uk*coeffs(:);


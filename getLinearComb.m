function [y] = getLinearComb(vecU,coeffs)
	%combines first k left singular vectors using coeffs and returns
	%coeffs must be k-length vector

	k=length(coeffs);
	Uk=reshape(vecU,[],k);
	y=Uk*coeffs(:);

	y=y(:);
	y=[y; zeros(length(vecU)-length(y),1)]

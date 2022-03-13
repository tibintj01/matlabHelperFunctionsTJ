function [isOnNormalSide]=pointIsOnNormalSide(p,n,q)
%p is point on plane, n is normal defining plane, q is query point
%returns 1 if q is on same side as normal of the plane, 0 otherwise
	numQpoints=size(q,1);
	isOnNormalSide=NaN(numQpoints,1);

	for qi=1:numQpoints
		isOnNormalSide(qi)=dot(q(qi,:)-p,n)>=0;
	end

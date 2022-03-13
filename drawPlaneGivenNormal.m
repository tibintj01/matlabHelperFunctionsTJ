function [sH] = drawPlaneGivenNormal(p,v,domainLimsX,domainLimsY,showNormal)

if(~exist('domainLimsX','var'))
	domainLimsX=[-1 1];
	domainLimsY=[-1 1];
end

if(~exist('showNormal','var'))
    showNormal=0;
end
x1=p(1);
y1=p(2);
z1=p(3);

domainX=linspace(domainLimsX(1),domainLimsX(2),10);
domainY=linspace(domainLimsY(1),domainLimsY(2),10);

if(sum(isnan(v))==length(v))
    sH=0;
    return
end

    w = null(v); % Find two orthonormal vectors which are orthogonal to v
   [P,Q] = meshgrid(domainX,domainY); % Provide a gridwork (you choose the size)
   X = x1+w(1,1)*P+w(1,2)*Q; % Compute the corresponding cartesian coordinates
   Y = y1+w(2,1)*P+w(2,2)*Q; %   using the two vectors in w
   Z = z1+w(3,1)*P+w(3,2)*Q;
   sH=surf(X,Y,Z);
   hold on
   alpha(0.2)
   if(showNormal)
   plot3DvectorFromTo(p,v,'k')
   end

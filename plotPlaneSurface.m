function [outputArg1,outputArg2] = plotPlaneSurface(planeNormalVec,pointOnPlane)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
N=100;

x1=pointOnPlane(1);
y1=pointOnPlane(2);
z1=pointOnPlane(3);


w = null(planeNormalVec(:)'); % Find two orthonormal vectors which are orthogonal to v
   [P,Q] = meshgrid(linspace(0,1,N)); % Provide a gridwork (you choose the size)
   X = x1+w(1,1)*P+w(1,2)*Q; % Compute the corresponding cartesian coordinates
   Y = y1+w(2,1)*P+w(2,2)*Q; %   using the two vectors in w
   Z = z1+w(3,1)*P+w(3,2)*Q;
   
   
   surf(X,Y,Z)



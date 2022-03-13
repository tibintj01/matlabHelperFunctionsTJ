function [xr] = rotate2Dvector(x,theta)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
%theta in radians

A=[cos(theta) -sin(theta); sin(theta) cos(theta)];

xr=A*x(:);

end


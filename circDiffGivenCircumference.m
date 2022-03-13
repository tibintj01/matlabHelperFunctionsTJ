function [circDiffVal] = circDiffGivenCircumference(val2,val1,circum)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


val1Radians=(val1/circum)*2*pi;
val2Radians=(val2/circum)*2*pi;

%angdiff on interval [-pi, pi] --> convert to [0,2*pi]
circDiffRad=angdiff(val1Radians,val2Radians);

circDiffVal=circDiffRad*circum/(2*pi);

end


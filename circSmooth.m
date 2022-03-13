function [circSmoothedVals] = circSmooth(values,windParam)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

wrappedValues=[values(:); values(:); values(:)];

wrappedSmoothedValues=smooth(wrappedValues,windParam);


startIdx=length(values)+1;
endIdx=startIdx+length(values)-1;

circSmoothedVals=wrappedSmoothedValues(startIdx:endIdx);




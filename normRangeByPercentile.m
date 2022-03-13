function [normRangeValues] = normRangeByPercentile(values,minPrc,maxPrc)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    minVal=prctile(values,minPrc);
    maxVal=prctile(values,maxPrc);
    
    normRangeValues=values-minVal;
    normRangeValues=normRangeValues/maxVal;
    
    
    
    


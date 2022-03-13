function [smoothed] = nanGaussSmoothRate(img)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    filtWidth = 3*2*3;
    filtSigma = 1.5*2*2*1.2*3;
    %filtSigma = 1;
    imageFilter=fspecial('gaussian',filtWidth,filtSigma);
    smoothed = nanconv(img,imageFilter, 'nanout','edge');
end


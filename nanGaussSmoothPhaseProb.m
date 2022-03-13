function [smoothed] = nanGaussSmoothPhaseProb(img)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    %filtWidth = 3*2;
    %filtSigma = 1.5*2;
    %filtSigma=0.6;
      %filtSigma=0.8*1.5;
         filtSigma=0.7;
      %filtSigma=1.2;
      
    filtWidth=2*ceil(2*filtSigma)+1;
    %filtSigma = 1;
    imageFilter=fspecial('gaussian',filtWidth,filtSigma);
    smoothed = nanconv(img,imageFilter, 'nanout','edge');
      %smoothed = nanconv(img,imageFilter, 'nanout');
end


function [pH] = plotCircXDistr(phaseBins,pDist,colorStr)
%UNTITLED2 Summary of this function goes here

pH=plot([ phaseBins(:); phaseBins(:)+360],[pDist(:); pDist(:)] ,colorStr,'LineWidth',3)
hold on





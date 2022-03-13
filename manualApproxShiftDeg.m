function [spikePhasesShifted] = manualApproxShiftDeg(spikePhases,fileBaseName)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    if(contains(fileBaseName,'atsby'))
             %bestShiftDeg=180;
                bestShiftDeg=100;
        elseif(contains(fileBaseName,'chilles'))
            %bestShiftDeg=-45; %45
             bestShiftDeg=-200; %45
            
        elseif(contains(fileBaseName,'uddy'))
            %bestShiftDeg=-45;
            bestShiftDeg=120;
         elseif(contains(fileBaseName,'icero'))
            bestShiftDeg=-230; %45
            
        else
             
             %bestShiftDeg=180;
              bestShiftDeg=0;
              %close all
              %continue
         end
        %first coarse then fine circ shift adjustment
        spikePhasesShifted=mod(spikePhases-bestShiftDeg,360);


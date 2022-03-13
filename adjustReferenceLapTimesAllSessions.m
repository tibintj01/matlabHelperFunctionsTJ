close all; clear all; clc
dataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles';
showPlots=0;
setTightSubplots_SpaceTime
%sessionNames={'Achilles_10252013','Achilles_11012013','Buddy_06272013', 'Cicero_09102014', 'Cicero_09172014','Gatsby_08282013'};

%Cicero_09172014 has only 1 strongly precessing cell used
sessionNames={'Achilles_10252013','Achilles_11012013','Buddy_06272013', 'Cicero_09102014','Gatsby_08282013'};

sessionNames={'Achilles_11012013','Buddy_06272013', 'Cicero_09102014','Gatsby_08282013'};

minCenterSeparationM=0.1; %counting similar place fields is like double counting, dilutes effect

for si=1:length(sessionNames)
    %try
    currSesName=sessionNames{si};
    
    [refUnitDataLeftward,refUnitDataRightward] = getRefUnitData(currSesName);
    
    for di=1:2
        if(di==1)
           currDirStr='rightward';
           refUnitData=refUnitDataRightward;
        else
           currDirStr='leftward';
           refUnitData=refUnitDataLeftward;
        end
    
        if(~isstruct(refUnitData) && isnan(refUnitData))
            continue
        end

        figure; plot(refUnitData.positionTimeAxis,refUnitData.directionSpecificStats.(currDirStr).lapNumberPerTime,'LineWidth',7)
        hold on
        
        lapStartTimes=refUnitData.lapStartTimesPerDir.(currDirStr);
        lapStopTimes=refUnitData.lapStopTimesPerDir.(currDirStr);
        
        numLaps=length(lapStartTimes);
        for li=1:numLaps
            plot([lapStartTimes(li) lapStartTimes(li)],ylim,'g-','LineWidth',3)
            plot([lapStopTimes(li) lapStopTimes(li)],ylim,'r-','LineWidth',3)
        end
        
        
        yyaxis right; plot(refUnitData.positionTimeAxis,refUnitData.positionPerTimeStep,'k')
    
    end


end

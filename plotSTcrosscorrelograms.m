close all; clear all; clc
dataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles';
showPlots=0;
setTightSubplots_SpaceTime
%sessionNames={'Achilles_10252013','Achilles_11012013','Buddy_06272013', 'Cicero_09102014', 'Cicero_09172014','Gatsby_08282013'};

%Cicero_09172014 has only 1 strongly precessing cell used
sessionNames={'Achilles_10252013','Achilles_11012013','Buddy_06272013', 'Cicero_09102014','Gatsby_08282013'};

halfLineHeight=0.01; %raster line height

for si=1:length(sessionNames)
    %try
    currSesName=sessionNames{si};
    trackFilePath=sprintf('%sTrackProps.mat',currSesName);
    currTrackData=load(trackFilePath);
    positionTimeAxis=currTrackData.positionTimeAxis;
    
    rowColInSeqMatrix=currTrackData.rowColInSequenceMatrix;
    
    %currTrackData=currTrackData.trackPopInfo

    for di=1:2
        if(di==1)
           currDirStr='rightward';
        else
           currDirStr='leftward';
        end
        if(~isfield(rowColInSeqMatrix,currDirStr))
            continue
        end
        currOverlapMatrix=rowColInSeqMatrix.(currDirStr);
        unitNames=fieldnames(currTrackData.trackPopInfo.(currDirStr)); %RELYING ON STABLE NUMBERING OF FIELDNAMES
        
        for fi=1:size(currOverlapMatrix,1)
            for fii=1:size(currOverlapMatrix,2)
                
                if(currOverlapMatrix(fi,fii)==0)
                    continue
                end
                
                currUnitName=unitNames{fi};
                overlappingNextUnitName=unitNames{fii};
                
                currUnitStruct=currTrackData.trackPopInfo.(currDirStr).(currUnitName);
                nextUnitStruct=currTrackData.trackPopInfo.(currDirStr).(overlappingNextUnitName);
                saveFileBaseName=sprintf('spacetimeCrossCorrelogram%sDir%d_%sVs%s',currSesName,di,currUnitName,overlappingNextUnitName);
                
                if(exist([saveFileBaseName '.png'],'file'))
                    %continue
                end
                %try
                    fH=figure
                    subplot(2,2,[1 3])
                    [currXcorrInfo]=getXcorrInfo(currUnitStruct,nextUnitStruct,positionTimeAxis,fH);

                     title(removeUnderscores(sprintf('%s, %s, %s vs %s',currSesName,currDirStr,currUnitName,overlappingNextUnitName)))
                     axis square
                     subplot(2,2,2)  
                      imshow(currUnitStruct.imagePath)
                      subplot(2,2,4)  
                       imshow(nextUnitStruct.imagePath)

                        setFigFontTo(18)
                         maxFig
                         print('-r75',saveFileBaseName,'-dpng')
                        close all

                        eval(sprintf('%sVs%sXcorrInfo = currXcorrInfo;',currUnitName,overlappingNextUnitName))

                    save(trackFilePath,sprintf('%sVs%sXcorrInfo',currUnitName,overlappingNextUnitName),'-append')
                  %{  
                catch ME
                    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
                    disp('SKIPPING:')
                    disp(ME.message)
                    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

                end
                    %}
            end
        end
    end
    %{
      catch ME
                    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
                    disp('SKIPPING:')
                    disp(ME.message)
                    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

                end
    %}
end
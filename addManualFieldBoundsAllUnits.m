close all; clear all; clc

processedDataDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData';
markedBoundsSaveDir='./hc3Images/phasePrecessionBounds';

touchDir(markedBoundsSaveDir)

unitDataDir=fullfile(processedDataDir,'unitSpikeInfoStructs');

unitFilePaths=getFilePathsRegex(unitDataDir,'*mat');

totalFieldCount=0;
%setTightSubplots

for fi=1:length(unitFilePaths)
    resaveImgFlag=0;
    
    if(mod(fi,50)==0)
        disp(fi)
    end
    currUnitFilePath=unitFilePaths{fi};
    currUnitStruct=load(currUnitFilePath);
    currUnitData=currUnitStruct.unitInfo;
    
    currUnitFileName=getFileNameFromPath(currUnitFilePath);
    
    
    %{
    if(strcmp(currUnitFileName,'ec013.156_Unit604Info.mat') || strcmp(currUnitFileName,'ec013.554_Unit806Info.mat') ...
            || strcmp(currUnitFileName,'ec013.574_Unit809Info.mat') || strcmp(currUnitFileName,'ec013.599_Unit712Info.mat'))
        resaveImgFlag=1;
        
    end
    %}
    %{
    posData=load(currUnitData.positionInfoForSessionSavePath);
    trackLengthM=posData.trackLengthM;
    
    rightSpikePositions=currUnitData.rightSpikePositions;
    rightSpikePhases=currUnitData.rightSpikePhases;
    leftSpikePositions=currUnitData.leftSpikePositions;
    leftSpikePhases=currUnitData.leftSpikePhases;
    %}
    if(~isfield(currUnitStruct,'rightwardFieldStartEndM') || mod(length(currUnitStruct.rightwardFieldStartEndM),2)~=0)
        plotHC3PlaceFieldsAndThetaPhases(currUnitData,'betterRightward')

        [rightwardFieldStartEndM,~] = ginput
        
        assert(mod(length(rightwardFieldStartEndM),2)==0)
        
        save(currUnitFilePath,'rightwardFieldStartEndM', '-append')
    elseif(~isempty(currUnitStruct.rightwardFieldStartEndM))
        
        plotHC3PlaceFieldsAndThetaPhases(currUnitData,'betterRightward')
        
        subplot(2,1,2)
        hold on
        bounds=currUnitStruct.rightwardFieldStartEndM;
        
        for bi=1:length(bounds)
            if(mod(bi,2)==1)
                plot([bounds(bi) bounds(bi)], [0 360],'g--','LineWidth',5)
            else
                plot([bounds(bi) bounds(bi)], [0 360],'r--','LineWidth',5)
            end
        end
        
        saveas(gcf,fullfile(markedBoundsSaveDir,sprintf('%s_Cell%dRightwardPhasePrecessionBounds.tif',currUnitData.sessionName,currUnitData.unitIDnum)))
        
        totalFieldCount=totalFieldCount+length(currUnitStruct.rightwardFieldStartEndM);
    end
    
    if(~isfield(currUnitStruct,'leftwardFieldStartEndM') || mod(length(currUnitStruct.leftwardFieldStartEndM),2)~=0)
        plotHC3PlaceFieldsAndThetaPhases(currUnitData,'betterLeftward')

        [leftwardFieldStartEndM,~] = ginput
        
        assert(mod(length(leftwardFieldStartEndM),2)==0)
        
        save(currUnitFilePath,'leftwardFieldStartEndM', '-append')
    elseif(~isempty(currUnitStruct.leftwardFieldStartEndM))
        
        plotHC3PlaceFieldsAndThetaPhases(currUnitData,'betterLeftward')
        subplot(2,1,2)
        hold on
        bounds=currUnitStruct.leftwardFieldStartEndM;
        
        for bi=1:length(bounds)
            if(mod(bi,2)==1)
                plot([bounds(bi) bounds(bi)], [0 360],'g--','LineWidth',5)
            else
                plot([bounds(bi) bounds(bi)], [0 360],'r--','LineWidth',5)
            end
        end
        
        saveas(gcf,fullfile(markedBoundsSaveDir,sprintf('%s_Cell%dLeftwardPhasePrecessionBounds.tif',currUnitData.sessionName,currUnitData.unitIDnum)))
        
        totalFieldCount=totalFieldCount+length(currUnitStruct.leftwardFieldStartEndM);
    end
    
    %close all
end
totalFieldCount/2
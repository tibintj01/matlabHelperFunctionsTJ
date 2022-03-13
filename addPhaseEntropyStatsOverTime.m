close all; clear all; clc

processedDataDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData';
unitDataDir=fullfile(processedDataDir,'unitSpikeInfoStructs');

unitFilePaths=getFilePathsRegex(unitDataDir,'*mat');

%setTightSubplots
for fi=1:length(unitFilePaths)
    currUnitFilePath=unitFilePaths{fi};
    currUnitStruct=load(currUnitFilePath);
    currUnitData=currUnitStruct.unitInfo;
    
    posData=load(currUnitData.positionInfoForSessionSavePath);
    trackLengthM=posData.trackLengthM;
    
    rightSpikePositions=currUnitData.rightSpikePositions;
    rightSpikePhases=currUnitData.rightSpikePhases;
    leftSpikePositions=currUnitData.leftSpikePositions;
    leftSpikePhases=currUnitData.leftSpikePhases;
    
    %{
    [leftPhaseHPerPos,posBinCenters]=getPhaseDistEntropyByPos(leftSpikePositions, leftSpikePhases,trackLengthM );
    [rightPhaseHPerPos,posBinCenters]=getPhaseDistEntropyByPos(rightSpikePositions, rightSpikePhases,trackLengthM);
    %}
    
    
     [leftLocalMIPerPos,posBinCenters]=getMovingMutualInformation(leftSpikePositions, leftSpikePhases,trackLengthM );
    [rightLocalMIPerPos,posBinCenters]=getMovingMutualInformation(rightSpikePositions, rightSpikePhases,trackLengthM);
    
    
    Nshuffles=1;
     %Nshuffles=10000;
     %Nshuffles=5000;
     Nshuffles=500;
    shuffledLeftPhaseEntropies=NaN(length(posBinCenters),Nshuffles);
    shuffledRightPhaseEntropies=NaN(length(posBinCenters),Nshuffles);
    
    for si=1:Nshuffles
        %leftShufflePhases=leftSpikePhases(randperm(length(leftSpikePhases)));
        %rightShufflePhases=rightSpikePhases(randperm(length(rightSpikePhases)));
        
        leftShufflePhases=rand(size(leftSpikePhases))*360;
        rightShufflePhases=rand(size(rightSpikePhases))*360;
      
        
        %leftShufflePos=leftSpikePositions(randperm(length(leftSpikePositions)));
        %rightShufflePos=rightSpikePositions(randperm(length(rightSpikePositions)));
        
        
        
        %[leftPhaseHPerPos,~]=getPhaseDistEntropyByPos(leftSpikePositions, leftShufflePhases,trackLengthM );
        %[rightPhaseHPerPos,~]=getPhaseDistEntropyByPos(rightSpikePositions, rightShufflePhases,trackLengthM);
        [leftPhaseHPerPos,~]=getMovingMutualInformation(leftSpikePositions, leftShufflePhases,trackLengthM );
        [rightPhaseHPerPos,~]=getMovingMutualInformation(rightSpikePositions, rightShufflePhases,trackLengthM);
        
        
        
        shuffledLeftPhaseEntropies(:,si)= (leftPhaseHPerPos(:));
        shuffledRightPhaseEntropies(:,si)= (rightPhaseHPerPos(:));
        
        if(mod(si,100)==0)
            disp(si)
        end
    end
    
    leftSigThresh=prctile(shuffledLeftPhaseEntropies,95,2);
    rightSigThresh=prctile(shuffledRightPhaseEntropies,95,2);
    
    figure;
    subplot(2,1,1)
    yyaxis left
    plot(posBinCenters,rightPhaseHPerPos,'LineWidth',3)
    hold on
    %plot(xlim,[rightSigThresh rightSigThresh],'k--','LineWidth',3)
    plot(posBinCenters,rightSigThresh,'k--','LineWidth',3)
    
    ylim([0 0.1])
       %ylim([0 1])
    yyaxis right
    plot(rightSpikePositions,rightSpikePhases,'k.')
    ylim([0 360])
    
    subplot(2,1,2)
    yyaxis left
    plot(posBinCenters,leftPhaseHPerPos,'LineWidth',3)
    hold on
     %plot(xlim,[leftSigThresh leftSigThresh],'k--','LineWidth',3)
     plot(posBinCenters,leftSigThresh,'k--','LineWidth',3)
    ylim([0 0.1])
     %ylim([0 1])
    
    yyaxis right
    plot(leftSpikePositions,leftSpikePhases,'k.')
    ylim([0 360])
        
    %save(currUnitFilePath,'', '-append')
end

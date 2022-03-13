function [lfpFileNames,cellPropFileNames]= getFileNamesAroundCh(ptDir,sessionID,ch)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Description:
%  collect files corresponding to surrounding channels
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
closeChannels=getSurroundingChs(ch,ptDir);

dataDir='/Users/tibinjohn/sampleHumanData/MG49session3';
lfpDataDir='decimatedLFP-MatFiles/Fs2000Hz';
cellPropDir='cellProperties-MatFiles';

cellCount=0;
for i=1:length(closeChannels)
        chStr=getChStr(closeChannels(i));
        lfpChanDir=sprintf('concatChan%d',closeChannels(i));
        lfpFileNames{i}=getRegexFilePath(fullfile(dataDir,lfpDataDir,lfpChanDir),'lfpDecConcat.mat');
        cellPropFileNamesCh=getRegexFilePaths(fullfile(dataDir,cellPropDir),sprintf('%s*_cell_properties*mat',chStr));
        for j=1:length(cellPropFileNamesCh)
                cellCount=cellCount+1;
                cellPropFileNames{cellCount}=cellPropFileNamesCh{j};
        end
end
%lfpFileNames
%cellPropFileNames
%[chMap, chRowCols]=getChMap(ptDir);
%chRowCols(closeChannels)

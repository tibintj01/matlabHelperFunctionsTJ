function [] = createConcatLFPfromMatDir(lfpMatDir)
	%Tibin John, 9/13/17
	%Description: concatenates all lfp files in directory
	%Input: name of directory containing lfp mat files (extracted from ns5)
	%Output: Saves file containing all lfp for all lfp mat files in directory
	lfpFileInfo=dir(fullfile(lfpMatDir,'*lfp*original*.mat'));

% 20110615-094530-023lfp96_original.mat

	subject='MG49';


	omarCombineLFPChanIntoNS5(subject,channelID,numFiles,fileList,outputFilename,outputParentDir,deltaTBetweenFiles,processedDataDir);


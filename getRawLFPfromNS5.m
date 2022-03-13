function [concatLFP] = getRawLFPfromNS5(sessionIdx,ch)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Description:
%
%
%Input:
%
%
%Output:
%
%
%Author: Tibin John, tibintj@umich.edu
%Project directory name: /nfs/turbo/lsa-ojahmed/tibin/spikeDynamicsAnalysisTibin 
%Created on 2018-09-06
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%session_list_human_seizures_NEW_FOR_ELLEN_SORTED_DATA_2017_TJ
[dataInfo]=session_list_human_seizures_NEW_FOR_ELLEN_SORTED_DATA_2017fxn();

%human_seizure_directories

SEIZURE_DATA_HOME_DIR='/nfs/turbo/lsa-ojahmed/tibin/FOR_TIBIN_otherPts/seizures'
SEIZURE_NEX_DATA_HOME_DIR='/nfs/turbo/lsa-ojahmed/tibin/spikeDynamicsAnalysisTibin/ellenMG49nexWithNoise'

CLASSIFICATION_FILES_DIR='/nfs/turbo/lsa-ojahmed/tibin/tibinCellSzClassificationFiles'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%extract input variables from struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('computing getRawLFPfromNS5 fxn output.........')
if(~strcmp(dataInfo(sessionIdx).subject,'MG49_seizure43'))
	ns5FilesDir=fullfile(SEIZURE_DATA_HOME_DIR,dataInfo(sessionIdx).subject,'NS5');
else
	ns5FilesDir=fullfile(SEIZURE_DATA_HOME_DIR,'MG49-seizure43','NS5');
end

ns5FilePath=getRegexFilePath(ns5FilesDir,sprintf('*_ch%d.ns5',ch))

disp('reading in NS5')
tic
concatLFPData=openNSx(ns5FilePath,'read');
toc

%convert digital to uV  - see humanProcesLFP... script:
%"This value is stored in the NEVs (openNEVOutput.ElectrodesInfo.DigitalFactor)"
concatLFP=round(double(concatLFPData.Data)*0.249*10000)/10000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compute output variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%store output variables in struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


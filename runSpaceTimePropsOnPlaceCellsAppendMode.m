%dataPathInfo=load('exampleHPCCells/placeCellFilePaths.mat');
dataPathInfo=load('exampleHPCCells/placeCellFilePathsMinPeakRate10.mat');
dataPathsTable=dataPathInfo.placeCellPaths;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INTERPOLATE LFP THETA CYCLE TIMES BETTER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dataDir='/Users/tibinjohn/thetaSeq/code/spaceTimeProcessedDir/manuallySortedPlaceCellsAbove10Hz/strong';

dataFilePaths=getFilePathsRegex(dataDir,'*mat');

%saveDir='spaceTimeProcessedDir';
%touchDir(saveDir)
saveDir=dataDir;

testMode=0;

for i=1:length(dataFilePaths)
   %if(mod(i,100)==0)
	 i/size(dataFilePaths,1)
    %end
    %try
	%{
	if(testMode)
		fileName='LEM3216_S20190726184722.mat'
		cellNum=3;
	else
        	fileName=dataPathsTable{i,1}{1};
        	cellNum=dataPathsTable{i,2};
	end
	%}

	currFilePath=dataFilePaths{i};
	currFileName=getFileNameFromPath(currFilePath);
	
	fileRootName=currFileName(1:22);
	fileName=[fileRootName '.mat'];
	cellNum=str2num(getSubstrBtwnSubstrs(currFileName,'cellNum','_dir'));
    
    di=str2num(getSubstrBtwnSubstrs(currFileName,'_dir','_space'));
    


	

	if(~contains(fileName,'LEM3116_S2018080610414') &&  ~contains(fileName,'LEM3116_S2018081414152') && ~contains(fileName,'LEM3206_S2019080717363') ...
		&& ~contains(fileName,'LEM3216_S2019072618472') && ~contains(fileName,'LEM3216_S2019072618472') && ~contains(fileName,'LEM3216_S2019081217560') ...
		&& ~contains(fileName,'LEM3216_S2019081313525') && ~contains(fileName,'LEM3216_S2019082313561') && ~contains(fileName,'LEM3246_S2019071711582') ...
		&& ~contains(fileName,'LEM3246_S2019072511441'))
		%fileName
		%continue
	end


        [spikePerCycleInfo] = getPerCycleSpikeInfoClarkFile(fileName,cellNum);
        
        %for di=1:2
           saveFileName=fullfile(saveDir,sprintf('%s_cellNum%d_dir%d_spaceTimePhaseCycleInfo.mat',fileRootName,cellNum,di));
               if(exist(saveFileName,'file'))
            %data=load(saveFileName);
            %spikePerCycleInfo=data.spikePerCycleInfo;
               end 
                if(di==1)
                    dirFlag=-1;
                else
                    dirFlag=1;
                end

                [spaceTimePhaseInfo] = getSpaceTimePhaseInfo(spikePerCycleInfo,dirFlag);
                %save(fullfile(saveDir,sprintf('%s_cellNum%d_dir%d_spaceTimePhaseCycleInfo.mat',fileRootName,cellNum,di))...
                %save(saveFileName...
            %	 ,'spaceTimePhaseInfo','spikePerCycleInfo');
            %try
            if(strcmp(currFilePath,saveFileName)==1)
                save(saveFileName...
                ,'spaceTimePhaseInfo','spikePerCycleInfo','-append'); %to update instead of erase all
            else
                disp('')
                
            end
            %catch
            %       save(saveFileName...
            %    ,'spaceTimePhaseInfo','spikePerCycleInfo'); %make new

            %end

            %checkIfGoodPlaceCell
            %skip=0;
            %if(skip)
            %    continue
            %end

            %testPlot

            %saveas(gcf,fullfile(saveDir,sprintf('%s_cellNum%d_dir%d.tif',fileRootName,cellNum,di)));
        %end
   % catch ME
   %     disp(ME.message)
   % end
	if(testMode)
		break
	end
end

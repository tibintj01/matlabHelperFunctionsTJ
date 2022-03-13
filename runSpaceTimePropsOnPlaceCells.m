%dataPathInfo=load('exampleHPCCells/placeCellFilePaths.mat');
dataPathInfo=load('exampleHPCCells/placeCellFilePathsMinPeakRate10.mat');
dataPathsTable=dataPathInfo.placeCellPaths;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INTERPOLATE LFP THETA CYCLE TIMES BETTER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
appendMode=1;

dataDir='/Users/tibinjohn/thetaSeq/code/spaceTimeProcessedDir/manuallySortedPlaceCellsAbove10Hz/strong';

saveDir='spaceTimeProcessedDir';
touchDir(saveDir)

testMode=1;
testMode=0;

for i=1:size(dataPathsTable,1)
   %if(mod(i,100)==0)
	 i/size(dataPathsTable,1)
    %end
    %try
	if(testMode)
		fileName='LEM3216_S20190726184722.mat'
		cellNum=3;
	else
        	fileName=dataPathsTable{i,1}{1};
        	cellNum=dataPathsTable{i,2};
	end

	%if(~contains(fileName,'LEM3116_S201808'))
	%	continue
	%end

	

	if(~contains(fileName,'LEM3116_S2018080610414') &&  ~contains(fileName,'LEM3116_S2018081414152') && ~contains(fileName,'LEM3206_S2019080717363') ...
		&& ~contains(fileName,'LEM3216_S2019072618472') && ~contains(fileName,'LEM3216_S2019072618472') && ~contains(fileName,'LEM3216_S2019081217560') ...
		&& ~contains(fileName,'LEM3216_S2019081313525') && ~contains(fileName,'LEM3216_S2019082313561') && ~contains(fileName,'LEM3246_S2019071711582') ...
		&& ~contains(fileName,'LEM3246_S2019072511441'))
		%fileName
		%continue
	end

        fileRootName=fileName(1:(end-1));

        [spikePerCycleInfo] = getPerCycleSpikeInfoClarkFile(fileName,cellNum);
        
        for di=1:2
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
            try
                save(saveFileName...
                ,'spaceTimePhaseInfo','spikePerCycleInfo','-append'); %to update instead of erase all
            catch
                   save(saveFileName...
                ,'spaceTimePhaseInfo','spikePerCycleInfo'); %make new

            end

            %checkIfGoodPlaceCell
            %skip=0;
            %if(skip)
            %    continue
            %end

            %testPlot

            %saveas(gcf,fullfile(saveDir,sprintf('%s_cellNum%d_dir%d.tif',fileRootName,cellNum,di)));
        end
   % catch ME
   %     disp(ME.message)
   % end
	if(testMode)
		break
	end
end

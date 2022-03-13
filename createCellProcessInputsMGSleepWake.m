
%cellProp=getCellWaveformPropsOmar_Tibin('MG49_seizure36','MG49_seizure36_20110612-151513-010_011_012',1,1,'/nfs/turbo/lsa-ojahmed/FOR_TIBIN_otherPts/seizures/MG49_seizure36/NEX',pwd);


%dirNames={'MG49_seizure36','MG49-seizure43','MG49-seizure45','MG63_seizure1-4','BW9_seizure1-3','RIHE1'};

%descripNames={'MG49_seizure36_20110612-151513-010_011_012','MG49_seizure43_20110612-151513-024_025_026','MG49_seizure45_20110613-112805-068_069_070','MG63_seizures1-4','BW9_seizures1-3','RIHE1_06-094327-001_to_011'};

descripNames={'20110615-094530-023_024_025_026_030_031_032_plus_anesthesia','20090215-173256-002_003_004_005_006_007_009_013_015'};
dirNames={'MG49','MG29'};

%saveDir=sprintf('seizureCellWaveformProperties-%s',getTimeStamp());
saveDir=sprintf('sleepWakeWaveformProperties-%s',getTimeStamp());

%oldSaveDir='/home/tibintj/compiledDir_20-Nov-2017_19-50-56_getCellWaveformPropsOmar_Tibin/seizureCellWaveformPropertiesRaw-20-Nov-2017_19-07-32';

%6 inputs into function being parallelized
inputCellArray=cell(1,6);
%session_list_human_seizures_NEW_FOR_ELLEN_SORTED_DATA_2017_TJ
session_list_human_lfp_spike_relationship_for_students


cellNum=0;
for dirIdx=1:length(dirNames)
	currDirName=dirNames{dirIdx};
	currFileRootName=descripNames{dirIdx};
	%loadDirName=sprintf('/nfs/turbo/lsa-ojahmed/FOR_TIBIN_otherPts/seizures/%s/NEX',currDirName);
	loadDirName=sprintf('/nfs/turbo/lsa-ojahmed/FOR_TIBIN_otherPts/%s/NEX',currDirName);

	dirInfo=dir(fullfile(loadDirName,'*.nex'));

	nexFileNames={dirInfo.name};

	for nexFileIdx=1:length(nexFileNames)
		nexFilePath=fullfile(loadDirName,nexFileNames{nexFileIdx})
		%tic
		%disp('reading nex file.......')
		%nexInfo=readNexFile(nexFilePath);
		%disp('done')
		%toc
		
		%if(isfield(nexInfo,'neurons'))
		%	numCells=length(nexInfo.neurons);
		%else
		%	continue
		%end	
		sessIdx=getSessionIdxMG(nexFilePath);
		cellChannels=dataInfo(sessIdx).cellChannel;
		%chStr=getSubstrBtwnSubstrs(nexFileNames{nexFileIdx},'_ch','-0');
		chStr=getSubstrBtwnSubstrs(nexFileNames{nexFileIdx},'_ch','.nex');
		ch=str2num(chStr);

		if(isempty(ch))
			%if(sessIdx==5)
			%	fds
			%end
			continue
			
		end

		cellChIdxes=find(cellChannels==ch);

		maxNumCells=length(cellChIdxes);
		%for cellIdx=1:numCells
		for cellIdx=1:maxNumCells
			cellChIdx=cellChIdxes(cellIdx);	
			if(dataInfo(sessIdx).cellQuality(cellChIdx)<3)
				continue
			end
			currSuffix=convertNumberToLetter(cellIdx);
				if(exist('oldSaveDir','var'))		
					cellPropSaveFile = fullfile(oldSaveDir,sprintf('%g%s_cell_properties_%s.mat',ch, currSuffix,currDirName));
					if(exist(cellPropSaveFile,'file'))
						continue
					end
				end
				%if(ch==10 && cellIdx==1 && length(findstr('BW',currDirName)>0))
				%	fds
				%end
				%cellNum
				cellNum=cellNum+1;
				inputCellArray{cellNum,1}=currDirName;
				inputCellArray{cellNum,2}=descripNames{dirIdx};
				inputCellArray{cellNum,3}=ch;
				inputCellArray{cellNum,4}=cellIdx;
				inputCellArray{cellNum,5}=loadDirName;
				inputCellArray{cellNum,6}=saveDir;
		end
	end
end

inputCellArray
save('mgSleepWakeParallelInputCell.mat','inputCellArray')

function [] = concatenateLFPchanAcrossMatFiles(sessNum,dataDir,chan)

	disp('concatenating...')
	tic

	%dataDir should be directly above recording directories
	 session_list_human_lfp_spike_relationship_for_students

	 if(~exist('dataDir'))
		dataDir='C:/Users/tibin/spikeDynamicsAnalysisTibin/session3Mat';
	 end

	 if(~exist('sessNum'))
		sessNum=3;
	 end
	 numChannels=96; %general?
	subject=dataInfo(sessNum).subject;
	 recordingDirs=dataInfo(sessNum).ns5Filenames;

		 concatenatedLFP=[];
		 for fileNum=1:length(recordingDirs)
			recordingDir=recordingDirs{fileNum};
			 filePath=fullfile(dataDir,subject,recordingDir,sprintf('lfp%d_dec.mat',chan))
		 	 
			 if(exist(filePath,'file'))
		 	 	lfpData=matfile(filePath);
			 else
			 	disp('channel not found, skipping....')
				continue
			 end
			%lfp=round(lfpData.lfp*0.249*1000)/1000;%convert digital to uV  - see humanProcesLFP... script
			lfp=round(double(lfpData.lfp)*0.249*10000)/10000;%convert digital to uV  - see humanProcesLFP... script
			concatenatedLFP=[concatenatedLFP lfp(:)'];
		 end
		
		 if(~exist(fullfile(dataDir,sprintf('concatChan%d',chan)),'dir'))
		 	mkdir(fullfile(dataDir,sprintf('concatChan%d',chan)))
		 else
		 	delete(fullfile(dataDir,sprintf('concatChan%d',chan),'*.mat'));
		 end

		save((fullfile(dataDir,sprintf('concatChan%d',chan),'lfpDecConcat.mat')),'concatenatedLFP','-v7.3')

	toc
	disp('done')
	

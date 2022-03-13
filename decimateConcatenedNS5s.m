function [] = decimateConcatenedNS5s(concatDir,ns5FileName)
	%ns5FilePaths=getRegexFilePaths(concatDir,'*ns5');

	ns5RootName=ns5FileName(1:(end-4));
	ns5FilePath=fullfile(concatDir,ns5FileName);

	data=openNSx(ns5FilePath,'read','report');
	lfp=data.Data;
	
	decFactor=15;
	lfpDec = resample(double(lfp(:)),1, decFactor);

	save(fullfile(concatDir,sprintf('dec15x_%s.mat',ns5RootName)),'lfpDec')
	


		

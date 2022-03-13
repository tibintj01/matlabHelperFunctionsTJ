function [rawSegment] = getRawSegmentFromConcat(startTime,endTime,chanNum,concatDir)

	if(~exist('concatDir'))
		%disp('Using directory for MG49 session 3......')
		concatDir='C:\Users\tibin\spikeDynamicsAnalysisTibin\session3Mat';
	end

	filePath=fullfile(concatDir,sprintf('concatChan%d',chanNum), 'lfpOriginalConcat.mat');

	Fs=30000;

	startSample=round(startTime*Fs);
	endSample=round(endTime*Fs);

	lfpData=matfile(filePath);
	rawSegment=lfpData.concatenatedLFP(1,startSample:endSample);

	%rawSegment=lfp(startSample:endSample);






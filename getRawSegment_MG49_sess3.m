function [rawSegment] = getRawSegment_MG49_sess3(timeStart,timeEnd,chanNum)
	%7 files with 72035268 samples at 30000 Hz, 8th file with 56705240 samples at 30000Hz
	Fs=30000;
 	startSample=timeStart*Fs;
	endSample=timeEnd*Fs;
	numSamplesInFile=72035268; %~40 min

	startFilePos=floor(startSample/numSamplesInFile);

	parentDirName='C:/Users/tibin/spikeDynamicsAnalysisTibin/session3Mat/MG49';
	recordingFileNames={'20110615-094530-023','20110615-094530-024','20110615-094530-025','20110615-094530-026','20110615-094530-030','20110615-094530-031','20110615-094530-032','20110616-145122-001'}

	fullPath=fullfile(parentDirName,recordingFileNames{startFilePos+1},sprintf('lfp%d_original.mat',chanNum))

	rawLFPdata=load(fullPath);

	rawLFP=rawLFPdata.lfp;

	if(endSample <= length(rawLFP))
		rawSegment=rawLFP(startSample:endSample);
	else
		rawSegment=rawLFP(startSample:end);
	end
 



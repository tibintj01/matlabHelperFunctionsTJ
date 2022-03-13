function [ns5RootName,sampStart,sampEnd]=getFileAndSampleNum(chanNum,timeStart,timeEnd)
	%Tibin John, 9/20/17
	%Description: gets file name and sample numbers given channel and time interval
	metaTagDir='C:/Users/tibin/spikeDynamicsAnalysisTibin/session3Mat-concatenated/MG49/20110615-094530-023_024_025_026_030_031_032_plus_anesthesia';

	dirNames=strsplit(metaTagDir,'/');

	fileRootName=dirNames{end};
	sampleFileName=sprintf([fileRootName '_ch%d_metatags.mat'],chanNum);
	fileName=fullfile(metaTagDir,sampleFileName);


	%disp('Meta data file loading.....')
	m=load(fileName);
	metaTags=m.allMetaTags;

	Fs=metaTags{1}.SamplingFreq;

	fileLengths=zeros(size(metaTags));
	
	for f=1:length(fileLengths)
		fileLengths(f)=metaTags{f}.NumofPackets;
	end
	
	totalNumSamples=sum(fileLengths);

	%disp('filling file array.....')
	
	%takes 3 seconds....
	%sampleToFileIdx=zeros(totalNumSamples,1);

	%globalStartIdx=1;
	%for f=1:length(fileLengths)
	%	currFileEnd=globalStartIdx+fileLengths(f)-1;
	%	sampleToFileIdx(globalStartIdx:currFileEnd)=f;
	%	globalStartIdx=currFileEnd+1;
	%end

	%disp('getting sample numbers........')
	sampStartGlobal=timeStart*Fs;
	sampEndGlobal=timeEnd*Fs;

	startFileIdx=sampleToFileIdx(sampStartGlobal,fileLengths);
	endFileIdx=sampleToFileIdx(sampEndGlobal,fileLengths);

	if(startFileIdx~=endFileIdx)
		error('Error: Time range given spans files!')
	end
	
	ns5Filename=metaTags{startFileIdx}.Filename;
	ns5RootName=ns5Filename(1:(end-4));

	if(startFileIdx>1)
		prevFileIdx=startFileIdx-1;
		for subF=prevFileIdx:-1:1
			sampStartGlobal=sampStartGlobal-fileLengths(subF);
			sampEndGlobal=sampEndGlobal-fileLengths(subF);
		end
	end

	sampStart=sampStartGlobal;
	sampEnd=sampEndGlobal;

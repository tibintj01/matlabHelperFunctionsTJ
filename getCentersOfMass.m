function [leftCMBin,rightCMBin] = getCentersOfMass(cellPropMatPath,leftFloor,rightFloor,saveDir)
	%floor values should between 0 and 1 (fraction of peak/trough amplitude)
	tic
	rootName=getFileNameFromPath(cellPropMatPath);
	rootName=rootName(1:(end-4));

	cellProps=load(cellPropMatPath);

	waveforms=cellProps.waveforms;
	Fs=30000;

	szStart=getSeizureStartTimes(cellPropMatPath);
	szStartTime=szStart(1);
	goodSpikesB4sz=intersect(cellProps.goodSpikes,(find(cellProps.spikeTimes<szStartTime)));

	goodWaveforms=waveforms(:,goodSpikesB4sz);

	%size(goodWaveforms)

	%fds

	meanGoodWaveform=mean(waveforms(:,goodSpikesB4sz),2);

	numBins=size(meanGoodWaveform,1);
	splineDt=0.1;
	avgWavInterped=spline(1:numBins, meanGoodWaveform, 1:splineDt:numBins);
	goodWaveformsInterp=zeros(length(1:splineDt:numBins),size(goodWaveforms,2));
	leftCMBins=zeros(size(goodWaveforms,2));
	rightCMBins=zeros(size(goodWaveforms,2));
	leftCMVals=zeros(size(goodWaveforms,2));
	rightCMVals=zeros(size(goodWaveforms,2));

	for i=1:size(goodWaveforms,2)
		wavInterped = spline(1:numBins, goodWaveforms(:,i), 1:splineDt:numBins);
	 	wavInterpedSmooth = fastsmooth(wavInterped, 20, 3, 1);	
		goodWaveformsInterp(:,i)=wavInterpedSmooth;

		[minVal,minIdxInterp]=min(wavInterpedSmooth);
		[maxVal,maxIdxInterp]=max(wavInterpedSmooth);
		[leftDataIdxes,leftDataSegment]=getTopPercentAroundPeak(-wavInterpedSmooth,minIdxInterp,leftFloor);
		[rightDataIdxes,rightDataSegment]=getTopPercentAroundPeak(wavInterpedSmooth,maxIdxInterp,rightFloor);
		[leftCMIdx,leftCMVals(i)]=getCenterOfMass(leftDataSegment);
		[rightCMIdx,rightCMVals(i)]=getCenterOfMass(rightDataSegment);
		leftCMBins(i)=leftCMIdx+leftDataIdxes(1)-1;
		rightCMBins(i)=rightCMIdx+rightDataIdxes(1)-1;
		%convert back to trough
	end
	leftCMBin=nanmean(leftCMBins);
	rightCMBin=nanmean(rightCMBins);

	leftCMBin=leftCMBin(1);
	rightCMBin=rightCMBin(1);
	%plot(meanGoodWaveform)

	%[inflectVal,inflectBin]=getInflectionPoint(avgWavInterped);
	%zeroedInterpWav=avgWavInterped-inflectBin;
	%plot(zeroedInterpWav,'k')
	
	%[minVal,minIdxInterp]=min(avgWavInterped);
	%[maxVal,maxIdxInterp]=max(avgWavInterped);
	%[leftDataIdxes,leftDataSegment]=getTopPercentAroundPeak(-avgWavInterped,minIdxInterp,leftFloor);
	%[rightDataIdxes,rightDataSegment]=getTopPercentAroundPeak(avgWavInterped,maxIdxInterp,rightFloor);
	%[leftCMIdx,leftCMVal]=getCenterOfMass(leftDataSegment);
	%[rightCMIdx,rightCMVal]=getCenterOfMass(rightDataSegment);
	%leftCMBin=leftDataIdxes(leftCMIdx);
	%rightCMBin=rightDataIdxes(rightCMIdx);
	%convert back to trough
	%leftCMVal=-leftCMVal;
	

	if(exist('saveDir'))
		plot(avgWavInterped,'r')
		hold on
		plot(leftCMBins,leftCMVals,'ko')
		plot(rightCMBins,rightCMVals,'ko')
		%saveas(gcf, fullfile(saveDir,['AllGoodSpikes_CM.tif']))	
		saveas(gcf, [rootName '_CM_IndividualWaveforms.tif'])
	end
	toc	

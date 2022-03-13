function [fileIdx] = sampleToFileIdx(sampleNum,fileLengths)

	totalSamples=0;
	fileIdx=0;
	while(totalSamples<sampleNum)
		fileIdx=fileIdx+1;
		totalSamples=totalSamples+fileLengths(fileIdx);
	end

	

function [instSpecSlope]=getInstSlopeFromPowerSpec(powerSpec,freqs,startFreq,endFreq)

	startFreqIdx=getIdxOf(startFreq,freqs);
	endFreqIdx=getIdxOf(endFreq,freqs);

	x=freqs(startFreqIdx:endFreqIdx);
	y=powerSpec(startFreqIdx:endFreqIdx);
	    
	bm=robustfit(x,y);
	instSpecSlope=bm(2);
	%instSpecSlope=(powerSpec(endFreqIdx)-powerSpec(startFreqIdx))/(endFreq-startFreq);


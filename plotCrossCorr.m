function [] =plotCrossCorr(spikeTimes1,spikeTimes2,titleStr)
	tic
	%Fs=30000;
	%binarize at different sampling rate to capture synchrony...
	%at 1/Fs time scale
	Fs=2000;
	spikeTrain1=rasterize(spikeTimes1,Fs);
	spikeTrain2=rasterize(spikeTimes2,Fs);
	isAuto=0;

	if(length(spikeTrain1)==length(spikeTrain2) && max(spikeTrain1-spikeTrain2)==0)
		isAuto=1;
	end

	maxLagTime=0.05;

	disp('computing xcorr.....')
	maxLag=maxLagTime*Fs;
	%[xc,lags]=xcorr(spikeTrain1,spikeTrain2,maxLag);
	binWidth=1/Fs;
	lags=-maxLag:1:maxLag;
	xc=crosscorrelation(spikeTrain1,spikeTrain2,binWidth,lags);

	%take out 0-lag peak if auto-correlation
	if(isAuto)
		xc(xc==max(xc))=0;
	end

	figure
	bar(lags/Fs*1000,xc)
	
	xlabel('Spike train lag (msec)')
	ylabel('Overlapping spike count')
	if(isAuto)
        title(sprintf('%s ac, %.2f ms bin width',titleStr,(1/Fs)*1000))
    else
        title(sprintf('%s xc, %.2f ms bin width',titleStr,(1/Fs)*1000))
	end
	toc

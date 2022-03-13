function [inspSpecSlope]=getInstSpecSlopeFromLFP(lfp,Fs)
	disp('Calculating spectrogram............')
	window=1;
	%winstep=0.01;
	winstep=1;
	lfp=double(lfp);
	[spectrogram f t]=omarSpectrogram(lfp,Fs,window,winstep);
	
	spectrogram=spectrogram';

	figure('units','normalized','outerposition',[0 0 1 1])
	ax1=subplot(3,1,1)
	plot((1:length(lfp))/Fs,lfp)

	ax2=subplot(3,1,2)
	colormap(ax2,fire)
	startFreq=30;
	endFreq=50;
	imagesc(spectrogram,[min(spectrogram(endFreq,:)),max(spectrogram(startFreq,:))])
	
	tTicStep=round(length(t)/10);
	fTicStep=round(length(f)/50);
	set(gca, 'XTick', 1:tTicStep:length(t), 'XTickLabel', t(1:tTicStep:length(t)));
	set(gca, 'YTick', 1:fTicStep:length(f), 'YTickLabel', f(1:fTicStep:length(f)));
	inspSpecSlope=zeros(size(t));
	ylim([getIdxOf(0,f), getIdxOf(50,f)])


	disp('Getting LFP spectrum slope at each time..........')
	for timeIdx=1:10:length(t)
		timeIdx/length(t)
		inspSpecSlope(timeIdx)=getInstSlopeFromPowerSpec(spectrogram(:,timeIdx),f,startFreq,endFreq);
	end
	ax3=subplot(3,1,3)

	timeIdxes=1:10:length(t);
	plot(t(timeIdxes),inspSpecSlope(timeIdxes))
	
	linkaxes([ax1, ax2, ax3],'x')
	xlim([0 500])

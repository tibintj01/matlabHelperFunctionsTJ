function [S_zeroed] = getSgramSVDforChNum(ch)

	[decLFP,decFs]=getDecLFPforCh(ch);

	disp('computing spectrogram....')
	tic
	%window=1;
	%winstep=0.05;
	window=30;
	winstep=30;
	%window=10;
	%winstep=10;
	%window=2;
	%winstep=2;
	tapers=[1 1];
	fpass=[0.5 120];
	[S f t]  =omarSpectrogram(decLFP,decFs,window,winstep,tapers,fpass,0);

	size(f)
	%fds
	size(t)

	%originally in different t as different rows
	size(S)
	S=S';

	%fds

        %startFreq=30;
        %endFreq=50;
        %imagesc(log(S))

	plotAndSaveFigures=1;
	if(plotAndSaveFigures)
		omarPcolor(t,f,log(S))
		%omarPcolor(t,f,S)
		%,[min(S(endFreq,:)),max(S(startFreq,:))])
		colormap(parula)
		xlabel('Time (sec)')
		ylabel('LFP Frequency (Hz)')
		shading flat
		caxis([-4 7])
		%caxis([exp(-4) exp(7)])
		c=colorbar
		c.Label.String='ln(Power)';
		title(sprintf('MG49 LFP spectrogram, Ch %d',ch))
		saveas(gcf,sprintf('MG49-LFP-LogSpectrogram-Ch%d.tif',ch))


		%subtract mean column
		n=size(S,2);
		idxCut=10000/winstep;
		S1=S(:,1:idxCut);
		S2=S(:,idxCut:end);
		n1=size(S1,2);
		n2=size(S2,2);

		meanSpec1=S1*ones(n1,1)*ones(n1,1)'/n1;
		meanSpec2=S2*ones(n2,1)*ones(n2,1)'/n2;

		figure
		plot(f,log(meanSpec1))
		hold on
		plot(f,log(meanSpec2))
		legend('Before 2.7hrs','After 2.7hrs')
		xlabel('LFP Freq. (Hz)')
		ylabel('ln(Power)')
		title(sprintf('MG49 LFP spec avg, Ch %d',ch))
	
		saveas(gcf,sprintf('MG49-LFP-LogSpecAvg-Ch%d.tif',ch))
	end
	%subtract mean column
	n=size(S,2);
		meanSpec=S*ones(n,1)*ones(n,1)'/n;

	S_zeroed=S-meanSpec;


	toc

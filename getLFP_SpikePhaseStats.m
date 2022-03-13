function [spikePhaseMean,spikePhaseKappa]=getLFP_SpikePhaseStats(spikeTimes,lfp,lfpTime,lowFreq,highFreq,titleStr,subPHandle,tenByTenPosition,maxKappaColor)
%function [spikePhaseMean,spikePhaseKappa,spikePhase,phaseBinEdges]=getLFP_SpikePhaseStats(spikeTimes,lfp,lfpTime,lowFreq,highFreq,titleStr,subPHandle,tenByTenPosition,maxKappaColor)
		%end

	%%establish settings
	makePlots=1;
	makePlotsHilbert=0;
	%maxKappaColor=0.6;
	lfpTimeAxis=lfpTime;	
	lfpFs=1/(median(diff(lfpTime)));

	%adapted from Omar's phase analysis code 
	filterOrder = 2; % go with 2 most of the time
	zscoreFilteredLFP = 1; % whether or not to zscore the filtered LFP

	
	%%filter LFP in specified band and convert to phase
	filtLFP = filterLFP(lfp, lfpFs, lowFreq, highFreq, filterOrder, zscoreFilteredLFP);

	%get phase of each bin in band-filtered lfp
	hilbertLFP = hilbert(-filtLFP); % NOTE THE INVERSION TO GET 180=MINIMA
	phaseLFP = angle(hilbertLFP);


	%% Convert phase from radians to degrees
	% phase ranges from -pi to +pi... make it range from 0 to 360
	phaseLFP = ((phaseLFP/pi) + 1)/2 * 360;
	

	if(makePlotsHilbert)
		lfpInd = find(lfpTimeAxis >= 1000.5 & lfpTimeAxis < 1001.5);
		figure;
		hold on;
		plot(lfpTimeAxis(lfpInd), zscore(filtLFP(lfpInd)), 'r');
		plot(lfpTimeAxis(lfpInd), zscore(phaseLFP(lfpInd)), 'b');
		title(titleStr)
		saveas(gcf,['PhaseSample-' titleStr '.tif'])
	end

 	% convert from spikeTimes to lfpBins
	spikeLFPBins = round(spikeTimes * lfpFs);

	%get phase of each spike
	spikePhase = phaseLFP(spikeLFPBins);
	
	%%calculate phase statistics of spikes
	spikePhaseMean = mod(rad2ang(circ_mean(ang2rad(spikePhase(:)))),360);
	[spikePhaseRTestP spikePhaseRTestZ] = circ_rtest(ang2rad(spikePhase(:)));
	%resultant vector lengths
	spikePhaseR = circ_r(ang2rad(spikePhase(:)));
	spikePhaseKappa   = circ_r2kappa(spikePhaseR(:))	
	kappaClrMap = copper(80);
	kappaIdx=round((spikePhaseKappa)*80/maxKappaColor);
	kappaIdx=min(80,kappaIdx);
	
	if(makePlots)
		phaseBinEdges = 0:20:360;
		phaseBinCenters = [10:20:350];
		spikeCount = histc(spikePhase, phaseBinEdges);
		spikePercent = spikeCount(1:end-1) ./ sum(spikeCount) * 100;

		%hold on;
		%plot(phaseBinCenters, spikePercent, 'x-', 'linewidth', 1, 'color', [0.1 0.1 0.1]);
		%subplot(10,10,1)
		%figure(10);
		%ha = tight_subplot(10,10,[.01 .03],[.1 .01],[.01 .01]);
		%axes(ha(tenByTenPosition))
		figure(10)	
		
		try
			axes(subPHandle(tenByTenPosition))
			polarhistogram(ang2rad(spikePhase),ang2rad(phaseBinEdges),'FaceColor',kappaClrMap(kappaIdx,:));
			%[tout,rout]=rose(ang2rad(spikePhase),length(phaseBinEdges));
			%polar(tout, rout);
			%[xout, yout] = pol2cart(tout, rout);
			%set(gca, 'nextplot', 'add');
			%fill(xout, yout, kappaClrMap(kappaIdx,:));
			set(gca,'FontSize',3)
			%figure(10)
			%title(['PhaseHistogram-' titleStr ', 180=Trough'])
			%set(get(gca,'title'),'Units','normalized','Position',[0.5 -0.1])
			thetaticklabels({})
			drawnow
		catch
			disp('skipping plot.....')
		end
		%title(sprintf('%.1f',spikePhaseKappa))
		%polarhistogram(ang2rad(spikePhase),ang2rad(phaseBinEdges));
		%axes('Color','none','XColor','none','YColor','none')
		%xlabel('Phase (degrees)');
		%ylabel('% of Spikes');
		%xlim([0 360]);
	end

	%saveas(gcf,['PhaseHistogram-' titleStr '.tif'])
	saveDir='/nfs/turbo/lsa-ojahmed/classifiedMG49-sleepWake';

	print(gcf,fullfile(saveDir,['PhaseHistogram-' titleStr]),'-dtiff','-r600')
        print(gcf,fullfile(saveDir,['PhaseHistogram-' titleStr]),'-depsc','-r600')
        %print(gcf,['PhaseHistogram-' titleStr],'-dpdf','-r600')

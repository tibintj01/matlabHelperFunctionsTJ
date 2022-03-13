function [] = plotSpikeHistHeatMapNex(waveforms,Fs,titleStr,figH,includeCbar,cmapName,clims)
        %[spikeHistHeatmap,tBins,vBins]=getSpikeHistHeapMap(waveforms,-0.04,0.04);
        %[spikeHistHeatmap,tBins,vBins]=getSpikeHistHeapMap(waveforms,-0.08,0.08);
        [spikeHistHeatmap,tBins,vBins]=getSpikeHistHeapMap(waveforms,-0.06,0.06);

	blackBgrnd=1;
	spikeHistHeatmap(spikeHistHeatmap==0)=NaN;	

        tVals=tBins/Fs;

	if(blackBgrnd)
		set(gca, 'color', [0 0 0]);
		hold on; 	
	end

        omarPcolor(tVals*1e3,vBins*1e3,spikeHistHeatmap,figH)
	shading flat

	%{
        xlabel('Time (ms)')
	if(strcmp(titleStr,'1'))
       		 ylabel('\muV')
        end

	xlim([0 1.6])
	xticks([0 1.6])
	xticklabels({0,1.6})
        %colormap jet
	%}
        %colormap(fire)
        %colormap(parula)
        %colormap(jet)
        %colormap(winter)
        %cBar1=colorbar
	if(exist('clims'))
		caxis(clims)	
	else
		if(~isnan(prctile(spikeHistHeatmap(:),80)))
			caxis([0 prctile(spikeHistHeatmap(:),90)*1.5])
		end
	end

	if(includeCbar)
		northColorBar
		currCAxis=caxis;
		ylabel(cbar1,sprintf('%d  Count  %d',currCAxis(1),currCAxis(2)))
	       set(cbar1,'YTick',[])
	end

	if(~exist('cmapName'))
		colormap(gca,copper)
	else
		colormap(gca,eval(cmapName))
	end
	box off


	%if(~isnan(prctile(spikeHistHeatmap(:),80)))
	%	caxis([0 prctile(spikeHistHeatmap(:),80)])
	%end
	%if(~isnan(prctile(spikeHistHeatmap(:),90)))
	%	caxis([0 prctile(spikeHistHeatmap(:),90)])
	%end
	      
	title(titleStr)

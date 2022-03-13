function [] = plotSpikeHistHeatMap(waveforms,Fs,titleStr,figH)
        [spikeHistHeatmap,tBins,vBins]=getSpikeHistHeapMap(waveforms);
        tVals=tBins/Fs;
        omarPcolor(tVals*1e3,vBins,spikeHistHeatmap,figH)
        xlabel('Time (ms)')
        ylabel('V (uV)')
        shading flat
        colormap jet
        cBar1=colorbar
        ylabel(cBar1,'Count')
        title(titleStr)

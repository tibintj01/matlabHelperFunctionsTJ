function [meanVtoP,meanHalfTrough]=getMeanShapeParameters(waveforms,Fs)

		 cellAvgWave=nanmean(waveforms,1);
                meanVtoP=getValleyToPeakTime(cellAvgWave,Fs);
                %maxSlope=getMaxSlope(cellAvgWave,Fs);
                meanHalfTrough=getHalfTroughWidth(cellAvgWave,Fs);

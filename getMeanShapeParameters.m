function [meanVtoP,meanHalfTrough]=getMeanShapeParameters(waveforms,Fs)

                vToP=zeros(size(waveforms,1),1);
                halfTrough=zeros(size(waveforms,1),1);

		numWaves=size(waveforms,1);
                for w=1:numWaves
                        waveform=waveforms(w,:);
                        vToP(w)=getValleyToPeakTime(waveform,Fs);
                        %maxSlope=getMaxSlope(waveform,Fs);
                        %halfPeak(w)=getHalfPeakWidth(waveform,Fs);
                        halfTrough(w)=getHalfTroughWidth(waveform,Fs);
                end
                meanVtoP=mean(vToP);
                meanHalfTrough=mean(halfTrough);

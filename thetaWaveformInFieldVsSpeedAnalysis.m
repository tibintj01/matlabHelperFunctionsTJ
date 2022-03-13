reload=0;

if(reload)
    close all; clear all; clc;
     tic;
    spikeDataPerField=load('allUnitInFieldGoodSpeedSpikePhasesAndTimes.mat');
    toc;  
end

 totalFieldCount=spikeDataPerField.totalFieldCount;
 thetaWaveformLFPzPerSpikePerField=spikeDataPerField.thetaWaveformLFPzPerSpikePerField;
 spikeAvgSpeedInTraversalsPerField=spikeDataPerField.spikeAvgSpeedInTraversalsPerField;
 
 %speedBinEdges=0.1*(1:8);
  speedBinEdges=0.1*(1:7);
 speedBinCenters=edgesToBins(speedBinEdges);
 
 numSpeedBins=length(speedBinCenters);
 
 timeEdges=linspace(0,0.14,50);
 %lfpZEdges=linspace(-3,3,200);
  lfpZEdges=linspace(-2,2,100);
 
 lfpVsTimePerSpeedBin=cell(numSpeedBins,1);
 
 setTightSubplots
 
 fH=figure;
 for fi=1:totalFieldCount
     fi
     currFieldThetaWaveformLFPzPerSpike=thetaWaveformLFPzPerSpikePerField{fi};
     currFieldSpikeAvgSpeedInTraversals=spikeAvgSpeedInTraversalsPerField{fi};
     
     for si=1:length(currFieldSpikeAvgSpeedInTraversals)
         
        currSpikeTravSpeed=currFieldSpikeAvgSpeedInTraversals(si);
        currSpeedBinIdx=discretize(currSpikeTravSpeed,speedBinEdges);
        
        if(isnan(currSpeedBinIdx))
            continue
        end
        
        currSpikeThetaWaveformLFPz=currFieldThetaWaveformLFPzPerSpike{si};
        
        if(length(lastLFPadded(:))==length(currSpikeThetaWaveformLFPz(:)))
            sameWaveformAsLast=max(abs(lastLFPadded(:)-currSpikeThetaWaveformLFPz(:)))>0;
        else
            sameWaveformAsLast=0;
        end
       
        if(~exist('lastLFPadded','var') || ~sameWaveformAsLast)
            lfpVsTimePerSpeedBin{currSpeedBinIdx}=[lfpVsTimePerSpeedBin{currSpeedBinIdx};currSpikeThetaWaveformLFPz];
        %else
        %    disp('repeat waveform')
            
        end
        
        lastLFPadded=currSpikeThetaWaveformLFPz;
         
     end
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %plot
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     if(mod(fi,100)==0)
         for spi=1:numSpeedBins
             subplot(numSpeedBins,1,spi)
            %plot(currFieldLowSpeedSpikeTimeFracs,currFieldLowSpeedSpikePhases,'b.')
            %getJointDistr(x,y,xEdges,yEdges,fH,withSmooth,cmapName)
            if(~isempty(lfpVsTimePerSpeedBin{spi}))
                getJointDistr(lfpVsTimePerSpeedBin{spi}(:,1)*1000,lfpVsTimePerSpeedBin{spi}(:,2),timeEdges*1000,lfpZEdges,fH,0,jet)
            end
            axis square
            xlabel('Time in cycle (msec)')
            ylabel('LFP (Z)')
            title(sprintf('%.2f < traversal speed < %.2f',speedBinEdges(spi),speedBinEdges(spi+1)))
         end
        drawnow
        maxFig
        disp('')
        
     end
 end
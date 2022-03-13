function [] = makeLocalMEAmovieForCh(ptDir,sessionNum,ch) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Description: makeMEAmovie.m
%
%
%Preconditions:
%
%
%Effects:
%
%
%Author: Tibin John, tibintj@umich.edu
%Project directory name: /nfs/turbo/lsa-ojahmed/tibin/processedHumanData/MG49/sessionID-3 
%Created on 2018-06-24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fs=2000;
sampPeriod=1/Fs;	
chStr=getChStr(ch);
maxFreq=80;
%maxFreq=30;
close all
maxR=8;

v=VideoWriter(sprintf('SpikeRateAndLockingVsLocalLFPRhythmsMovie-%s-%d-Ch%s.avi',ptDir,sessionNum,chStr));
v.FrameRate=2;
%v.FrameRate=4;
open(v);
figure(1)
maxFig
numPPCclrs=80;
ppcClrMap=copper(numPPCclrs);
maxPPC=0.3;
lastNumSpikes=0;

phaseBinEdges = 0:20:360;
phaseBinCenters = [10:20:350];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Task #1
%Description:
%  collect files corresponding to surrounding channels
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[lfpFileNames,cellPropFileNames]=getFileNamesAroundCh(ptDir,sessionNum,ch);


%[chMap, chRowCols]=getChMap(ptDir);
%chRowCols(closeChannels)

%load lfp for current channel
%EDIT TO REFLECT DATA LOCATION
dataDir='/Users/tibinjohn/sampleHumanData/MG49session3';
lfpDataDir='decimatedLFP-MatFiles/Fs2000Hz';
lfpChanDir=sprintf('concatChan%d',ch);
centerLfpFileName=getRegexFilePath(fullfile(dataDir,lfpDataDir,lfpChanDir),'lfpDecConcat.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Task #2
%Description:
% get (compute or load) spectrograms for this cluster of channels
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
centerLfpData=load(centerLfpFileName);
centerLfp=centerLfpData.concatenatedLFP;
[tSgram fSgram Sgram] = getSpectrogram(centerLfp,Fs);
specgramMatrix=NaN(length(fSgram),length(tSgram),length(lfpFileNames));
%lfpTimeAxis=sampPeriod:sampPeriod:(sampPeriod*length(lfp));
lfpTimeAxis=sampPeriod:sampPeriod:(sampPeriod*length(centerLfp));

specDataFileName=sprintf('%s_SessionID%d_%s_LocalSpectrogramsMatrix.mat',ptDir,sessionNum,chStr);

if(~exist(specDataFileName))
    for chIdx=1:length(lfpFileNames)
        currLFPFileName=lfpFileNames{chIdx};
        lfp=load(currLFPFileName);
        [tSgram fSgram Sgram] = getSpectrogram(lfp,Fs);
        specgramMatrix(:,:,chIdx)=Sgram;
        save(specDataFileName,'specgramMatrix')
    end
else
    specgramMatrixData=load(specDataFileName);
    specgramMatrix=specgramMatrixData.specgramMatrix;
end
%figure
%omarPcolor(tSgram,fSgram,log(specgramMatrix(:,:,1)))
%shading flat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Task #3
%Description:
% get (compute or load) band-filtered time traces for this cluster of channels
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for chIdx=1:length(closeChannels)
%	ch=closeChannels(chIdx);

	%CHANGE THIS TO REFLECT DATA LOCATION
	saveDir=pwd;
	fpass=[0.5 4];
	[DeltaLFPprops] = getFiltAndPhaseLFP(centerLfp,ch,fpass,saveDir);
	
	%fpass=[5 10];
	fpass=[3.5 8.5];
	[ThetaLFPprops] = getFiltAndPhaseLFP(centerLfp,ch,fpass,saveDir);
	
	fpass=[9.5 14.5];
	[AlphaLFPprops] = getFiltAndPhaseLFP(centerLfp,ch,fpass,saveDir);
	
	%fpass=[30 50];
	%[GammaLFPprops] = getFiltAndPhaseLFP(ch,fpass,saveDir);

%end

lfp=load(lfpFileNames{1});

[tSgram fSgram Sgram] = getSpectrogram(lfp,Fs);
specgramMatrix=NaN(length(fSgram),length(tSgram),length(lfpFileNames));
lfpTimeAxis=sampPeriod:sampPeriod:(sampPeriod*length(lfp));

specDataFileName=sprintf('%s_SessionID%d_%s_LocalSpectrogramsMatrix.mat',ptDir,sessionNum,chStr);

if(~exist(specDataFileName))
    for chIdx=1:length(lfpFileNames)
        currLFPFileName=lfpFileNames{chIdx};
        lfp=load(currLFPFileName);
        [tSgram fSgram Sgram] = getSpectrogram(lfp,Fs);
        specgramMatrix(:,:,chIdx)=Sgram;
        save(specDataFileName,'specgramMatrix')
    end
else
    specgramMatrixData=load(specDataFileName);
    specgramMatrix=specgramMatrixData.specgramMatrix;
end
%figure
%omarPcolor(tSgram,fSgram,log(specgramMatrix(:,:,1)))
%shading flat



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Task #4
%Description:
%plot power spectra of closest channels as a function of time
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%disp('computing output.........')

figure(1)
		uberTitle(sprintf('Local Spike Times And Phases vs Local LFP Rhythms Movie, %s-Session%d, Ch%s',ptDir,sessionNum,chStr))
		%axBar=subplotCartesian(0.5,0.5,0.01,0.01);	
		axBar=subplotCartesian(0.35,0.5,0.01,0.01);	
		%axBar=subplotCartesian(0,0.4,0.01,0.01);	
		colormap(copper)
		%position is from left, from bottom, width,
		%cHandle=colorbar('Position',[0.95 0.03 0.007 0.9]) 
		%cHandle=colorbar('Position',[0.5 0.5 0.007 0.4])
		cHandle=colorbar('Position',[0.35 0.5 0.007 0.4])
		%cHandle=colorbar('Position',[0 0.4 0.007 0.4])
		%cHandle=colorbar('Position',[0.5 0.5 0.01 0.01])
		set(cHandle,'FontSize',8)
		caxis([0 maxPPC])
		 ylabel(cHandle,'PPC','FontSize',10)
		set(axBar, 'Visible','off')
%fds

halfSmoothWind=1
sSampPeriod=tSgram(2)-tSgram(1);
numSmoothWindHalf=round(halfSmoothWind/sSampPeriod);

numCells=length(cellPropFileNames);
for i=1:numCells
    cellPropsData=load(cellPropFileNames{i});
    cellProps{i}=cellPropsData;
end

%rasterWindTimeHalf=10;
%rasterWindTimeHalf=5;
rasterWindTimeHalf=2.5;
%phaseWindTimeHalf=30;
phaseWindTimeFull=30;
%phaseWindTimeFull=120;
%phaseWindTimeHalf=15;

frameCount=0;
lfpLength=length(centerLfp);
for tIdx=(numSmoothWindHalf+1):(length(tSgram)-numSmoothWindHalf-1)
    currTimeMin=tSgram(tIdx)-sSampPeriod/2;
    currTimeMax=tSgram(tIdx)+sSampPeriod/2;

	axisStartTime=max((currTimeMin-rasterWindTimeHalf),1/Fs);
	axisEndTime=min((currTimeMax+rasterWindTimeHalf),lfpLength/Fs);
	currTimeAxis=axisStartTime:(1/Fs):axisEndTime;
	currLFPindices=round(currTimeAxis*Fs);
	%currLFPindices=currLFPindices(currLFPindices>=1 & currLFPindices<=lfpLength);
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%plot local power spectra at each time step (with moving avg window) 
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %subplot(121)
	%subplotCartesian(0,0.5,1/2,0.5);
  	subplotCartesian(0,2/3,1/3,1/3);
    for chIdx=1:length(lfpFileNames)
        localPowerSpec=nanmean(specgramMatrix(:,(tIdx-numSmoothWindHalf):(tIdx+numSmoothWindHalf),chIdx),2);
        plot(fSgram,log(localPowerSpec))
        hold on
        ylim([-3 8])
        xlim([0 maxFreq])
    end
    title(sprintf('9 local spectra of Ch%d; Time=%.2f sec',ch, tSgram(tIdx)))
    xlabel('Freq (Hz)')
    ylabel('log(LFP Power)')
       hold off

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%plot center filtered LFPs in time window  
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	deltaAx=subplotCartesian(2/3,0,1/3,1/7);
	%axes(deltaAx)
	plot(currTimeAxis,DeltaLFPprops.filtLFP_Z(currLFPindices))
	hold on
	plot([tSgram(tIdx) tSgram(tIdx)],[-4 4],'k--')
	hold off
	xlim([-Inf Inf])
	ylim([-4 4])
	title('Center delta bandpass Z (0.5-4Hz)')	
	set(gca,'Xticklabel',[])
	
	thetaAx=subplotCartesian(2/3,1/5,1/3,1/7);
	%axes(thetaAx)
	plot(currTimeAxis,ThetaLFPprops.filtLFP_Z(currLFPindices))
	hold on
	plot([tSgram(tIdx) tSgram(tIdx)],[-4 4],'k--')
	hold off
	title('Center theta bandpass Z (3.5-8.5Hz)')	
	xlim([-Inf Inf])
	ylim([-4 4])
	set(gca,'Xticklabel',[])

		
	alphaAx=subplotCartesian(2/3,2/5,1/3,1/7);
	%axes(alphaAx)
	plot(currTimeAxis,AlphaLFPprops.filtLFP_Z(currLFPindices))
	hold on
	plot([tSgram(tIdx) tSgram(tIdx)],[-4 4],'k--')
	hold off
	title('Center alpha bandpass Z (9.5-14.5Hz)')
	xlim([-Inf Inf])
	ylim([-4 4])
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%plot local spike raster in time window  
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %subplot(1,2,2)
  subplotCartesian(2/3,2/3,1/3,1/3);
      plot([tSgram(tIdx) tSgram(tIdx)],[0 numCells],'k--')
	title(sprintf('Ch%d 9 local probe cells raster',ch));
	 %axes('Position',[1/2 1/5 1/2 1/5])
   spikeTimesInPhaseWindAllCells=[];

   for cellIdx=1:numCells
  	subplotCartesian(2/3,2/3,1/3,1/3);
        spikeTimes=cellProps{cellIdx}.spikeTimes;
        %spikeTimesInWind=spikeTimes(spikeTimes>=currTimeMin & spikeTimes<= currTimeMax);
        %spikeTimesInWind=spikeTimes(spikeTimes>=axisStartTime & spikeTimes<= axisEndTime);
        spikeTimesInWind=spikeTimes(spikeTimes>=axisStartTime & spikeTimes<=tSgram(tIdx));
        if(length(spikeTimesInWind)>0)
            if(cellProps{cellIdx}.isInterneuronCell==1)
                %plot([spikeTimesInWind spikeTimesInWind], [cellIdx-1 cellIdx],'r')
                plot(spikeTimesInWind,cellIdx,'ro','MarkerSize',3)
            elseif(cellProps{cellIdx}.isInterneuronCell==0)
                %plot([spikeTimesInWind spikeTimesInWind], [cellIdx-1 cellIdx],'b')
                 plot(spikeTimesInWind,cellIdx,'bo','MarkerSize',3)
            else
                %plot([spikeTimesInWind spikeTimesInWind], [cellIdx-1 cellIdx],'k')
                 plot(spikeTimesInWind,cellIdx,'ko','MarkerSize',3)
            end
        end
        hold on
	    xlabel('Time (sec)')
	    ylabel('Cell No.')
	    ylim([0 numCells])
	    %xlim([currTimeMin-rasterWindTimeHalf currTimeMax+rasterWindTimeHalf])
	    xlim([axisStartTime axisEndTime])


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%plot local spike phase histogram in time window  
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%spikeTimesInPhaseWind=spikeTimes(spikeTimes>=tSgram(tIdx)-phaseWindTimeHalf & spikeTimes<=tSgram(tIdx)+phaseWindTimeHalf);
	spikeTimesInPhaseWind=spikeTimes(spikeTimes>=(tSgram(tIdx)-phaseWindTimeFull) & spikeTimes<=tSgram(tIdx));
	spikeTimesInPhaseWindAllCells=[spikeTimesInPhaseWindAllCells;spikeTimesInPhaseWind(:)];

	%if(tSgram(tIdx)>100)
	%	fds
	%end


    end
	    hold off
	
	if(length(spikeTimesInPhaseWindAllCells)~=lastNumSpikes)
		lastNumSpikes=length(spikeTimesInPhaseWindAllCells);
		
		spikeTimeIndicesInPhaseWind=round(spikeTimesInPhaseWindAllCells*Fs);
		
		deltaAx2=subplotCartesian(0,0,1/3,1/7);
		spikeDeltaPhasesInWind=DeltaLFPprops.phaseLFP(spikeTimeIndicesInPhaseWind);
		ppcDelta=omarPPC(spikeDeltaPhasesInWind);
	
		%polarhistogram(ang2rad(mod(spikeDeltaPhasesInWind+90,360)),ang2rad(phaseBinEdges))		
		ppcIdx=max(1,round(ppcDelta*(numPPCclrs/maxPPC)));
		ppcIdx=min(ppcIdx,numPPCclrs);
		polarhistogram(ang2rad(mod(spikeDeltaPhasesInWind+90,360)),ang2rad(phaseBinEdges),'FaceColor',ppcClrMap(ppcIdx,:))		
		thetaticklabels({})
		rlim([0 maxR])
		rticks([maxR])
		title(sprintf('Delta spike phases (past %d sec)',phaseWindTimeFull))	

		thetaAx2=subplotCartesian(0,1/5,1/3,1/7);
		spikeThetaPhasesInWind=ThetaLFPprops.phaseLFP(spikeTimeIndicesInPhaseWind);
		%polarhistogram(ang2rad(spikeThetaPhasesInWind),ang2rad(phaseBinEdges))		
		%polarhistogram(ang2rad(mod(spikeThetaPhasesInWind+90,360)),ang2rad(phaseBinEdges))		
		ppcTheta=omarPPC(spikeThetaPhasesInWind);
		ppcIdx=max(1,round(ppcTheta*(numPPCclrs/maxPPC)));
		ppcIdx=min(ppcIdx,numPPCclrs);
		polarhistogram(ang2rad(mod(spikeThetaPhasesInWind+90,360)),ang2rad(phaseBinEdges),'FaceColor',ppcClrMap(ppcIdx,:))		
		thetaticklabels({})
		rlim([0 maxR])
		rticks([maxR])
		title(sprintf('Theta spike phases (past %d sec)',phaseWindTimeFull))	

		alphaAx2=subplotCartesian(0,2/5,1/3,1/7);
		spikeAlphaPhasesInWind=AlphaLFPprops.phaseLFP(spikeTimeIndicesInPhaseWind);
		%polarhistogram(ang2rad(mod(spikeAlphaPhasesInWind+90,360)),ang2rad(phaseBinEdges))		
		ppcAlpha=omarPPC(spikeAlphaPhasesInWind);
		ppcIdx=max(1,round(ppcAlpha*(numPPCclrs/maxPPC)));
		ppcIdx=min(ppcIdx,numPPCclrs);
		polarhistogram(ang2rad(mod(spikeAlphaPhasesInWind+90,360)),ang2rad(phaseBinEdges),'FaceColor',ppcClrMap(ppcIdx,:))		
		thetaticklabels({})
		rlim([0 maxR])
		rticks([maxR])
		title(sprintf('Alpha spike phases (past %d sec)',phaseWindTimeFull))	
	end	








    %pause(1)
    frameCount=frameCount+1;
    %localMEAmovie(frameCount) = getframe(gcf);
    currFrame = getframe(gcf);
    writeVideo(v,currFrame);
    %pause(0.1)
    %if(mod(tSgram(tIdx),rasterWindTimeHalf*2)<sSampPeriod)
    %    cla
    %end
    if(frameCount>500)
    %if(frameCount>20)
        %figure; 
	%maxFig
	%axes('Position',[0 0 1 1]); movie(localMEAmovie,1)
	close(v);
        fds
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Task #3
%Description:
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('computing output.........')



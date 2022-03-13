function [outputStruct] = saveDecLFPAndCellPropMatFromConcatNS5Ch(sessionID,ch)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Description:
%
%
%Input:
%
%
%Output:
%
%
%Author: Tibin John, tibintj@umich.edu
%Project directory name: /nfs/turbo/lsa-ojahmed/tibin/spikeDynamicsAnalysisTibin 
%Created on 2018-06-26
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
session_list_human_lfp_spike_relationship_for_students

disp(sprintf('loading/saving processed data for ch %d.......',ch))

originalFs=30000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%extract input variables from struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
concatSessionName=dataInfo(sessionID).descriptiveFilename;

concatSessionNameCell=concatSessionName;
if(strcmp('20090216-184530-025_026_027_REMOVED_anesthesia',concatSessionName))
	concatSessionNameCell='20090216-184530-025_026_027_plus_anesthesia';
end

ptDir=dataInfo(sessionID).subject;

%nexFileDir=sprintf('/nfs/turbo/lsa-ojahmed/tibin/FOR_TIBIN_otherPts/sorted/concatenated/%s/%s/chan%d',ptDir,concatSessionName,ch);

nexFileDir=sprintf('/nfs/turbo/lsa-ojahmed/tibin/FOR_TIBIN_otherPts/sorted/concatenated/%s/%s/chan%d',ptDir,concatSessionNameCell,ch);

if(strcmp(ptDir,'RIHE1'))
	nexFileDir='/nfs/turbo/lsa-ojahmed/tibin/FOR_TIBIN_otherPts/seizures/RIHE1/NEX';
end


nexFilePath=getRegexFilePath(nexFileDir,sprintf('*_ch%d.nex',ch));
if(strcmp(ptDir,'RIHE1'))
	nexFilePath=getRegexFilePath(nexFileDir,sprintf('*_ch%d-01.nex',ch));
end

sortedCellSpikeTimesData=nexReadFile(nexFilePath);

%RMS, MAD; safe distance from noise; is it getting wider
%document everything method; here's all the spikes that were left
%manual sorting
%quantify how many wrong shape
%template matching
%spikes that were extras 
%convinces 
%MG49 ch43,ch1 - good cell - rated 1 instead of 5
%first half of seizures highest firing rates
%missing because got wider - how many
%post-ictal is depol block - well accepted
%caveat is proper sorting
%oscillate around a membrane - minimum for any given time
%recovery from depol block is only scenario where amp goes up
%can throw in supplements
%if in the noise, can say don't know; question is do we _know_ if anything else missed
%core vs penumbra
%if every FS is still spiking - did it change significantly
%chlloride became excitatory
%spike sorting specific to seizure time; could be knocked out by ellipse
%cell prop feature space - just looking at carefully, 
%how to cluster when everyone is changing? Spike shape matches pre seizure
%If RS spike is like unsorted spike, and FS is narrower 
%putative unknown cell if not like any others; two versions of spike rate; putatively belonging
%wanted  - people want to say my model reproduces data features
%one spike - people say hypersynchronous why you can't sort vs depol
%can easily detect each spike and wave with simple threshold - see if FS and RS have different phases
%depol block,synaptic depression, chloride reversal  - what is most likely first and what conditions
%simulatanous off but not on

%summarizing across brian states, irrespesctive of brain region; temporal and rate, and space and shape
%don't need task 
%add in ecog; pick closest 
%spikes with ecog; delta might be interesting across whole array
%phase locking to ecog; ewell in anesthesia
%angle matters for making it accessible - if svd, method might be better; 
%delta traveling waves and alpha supression in human brain
egSpikeTimes=sortedCellSpikeTimesData.neurons{1}.timestamps;
size(egSpikeTimes)

ns5FilesDir=sprintf('/nfs/turbo/lsa-ojahmed/tibin/processedHumanData-concatenated/%s/%s',ptDir,concatSessionName);

if(strcmp(ptDir,'RIHE1'))
	ns5FilesDir='/nfs/turbo/lsa-ojahmed/tibin/FOR_TIBIN_otherPts/seizures/RIHE1/NS5';
end

ns5FilePath=getRegexFilePath(ns5FilesDir,sprintf('*_ch%d.ns5',ch))
metaTagsPath=getRegexFilePath(ns5FilesDir,sprintf('*_ch%d_metatags.mat',ch))


cellPropSaveDir=sprintf('/nfs/turbo/lsa-ojahmed/tibin/processedHumanData/%s/sessionID-%d/cellProperties-MatFiles',ptDir,sessionID);

touchDir(cellPropSaveDir)

metaTags=load(metaTagsPath)

disp('reading in NS5')
tic
concatLFPData=openNSx(ns5FilePath,'read');
toc

%convert digital to uV  - see humanProcesLFP... script:
%"This value is stored in the NEVs (openNEVOutput.ElectrodesInfo.DigitalFactor)"
concatLFP=round(double(concatLFPData.Data)*0.249*10000)/10000;

decFactor=15;
decFs=originalFs/decFactor;

decLFP=concatLFP(1:decFactor:end);

saveLFPDir=sprintf('/nfs/turbo/lsa-ojahmed/tibin/processedHumanData/%s/sessionID-%d/decimatedLFP-MatFiles/Fs2000Hz',ptDir,sessionID);

touchDir(saveLFPDir)

chStr=getChStr(ch);
save(fullfile(saveLFPDir,sprintf('DecLFP_%s_ch%s.mat',concatSessionName,chStr)),'decLFP','decFs')



%testing if spike times correspond to lfp visually across concatenations - yes they do!
%{
figure
for i=length(egSpikeTimes):-1:(length(egSpikeTimes)-20)
%for i=1:10
	%subplot(10,1,i)
	spikeLFPIdx=(egSpikeTimes(i))*originalFs;
	%plot(concatLFP((spikeLFPIdx-10):(spikeLFPIdx+32)))
	plot(concatLFP((spikeLFPIdx-100):(spikeLFPIdx+200)))
	hold on
end
fds
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%start cell waveform processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('computing output.........')
cellChannels=dataInfo(sessionID).cellChannel;

cellChIdxes=find(cellChannels==ch);

%fds

maxNumCells=length(cellChIdxes);
for cellIdx=1:maxNumCells
	cellChIdx=cellChIdxes(cellIdx);
		%if(dataInfo(sessIdx).cellQuality(cellChIdx)<3)
		%	continue
		%end
	disp(sprintf('getting cell %d out of %d waveform properties.....',cellIdx,maxNumCells))	
	%getCellWaveformPropsOmar_Tibin_PtSleepWake(ptDir,sessionID,concatSessionName, ch, cellIdx, nexFileDir,ns5FilesDir, cellPropSaveDir)
	getCellWaveformPropsOmar_Tibin_PtSleepWake(ptDir,sessionID,concatSessionNameCell, ch, cellIdx, nexFileDir,ns5FilesDir, cellPropSaveDir)

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%store output variables in struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


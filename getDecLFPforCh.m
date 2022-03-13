function [decLFP,decFs] = getDecLFPforCh(subj,chanNum,sessionID,useClosest)

	if(useClosest)
		chStrSame=getChStr(chanNum);
		localKappaFile=sprintf('/nfs/turbo/lsa-ojahmed/tibin/processedHumanData/%s/sessionID-%d/spatialProperties/%s-%d-Ch%sNeighborLFPLFPKappas.mat',subj,sessionID,subj,sessionID,chStrSame);
		localKappaData=load(localKappaFile);
		localKappaMatrix=localKappaData.localKappaMatrix;

		%deltaThetaAlphaMeanKappa=mean(localKappaMatrix(:,1:3),2);
		allBandMeanKappa=nanmean(localKappaMatrix(:,1:end),2);
		sameChanIdx=find(localKappaData.neighborChs==chanNum);
		%deltaThetaAlphaMeanKappa(sameChanIdx)=NaN;
		allBandMeanKappa(sameChanIdx)=NaN;
		%[maxKappaAvg,bestLocalChanIdx]=max(deltaThetaAlphaMeanKappa);
		[maxKappaAvg,bestLocalChanIdx]=max(allBandMeanKappa);
		bestLocalChanNum=localKappaData.neighborChs(bestLocalChanIdx)
	
		%localKappaData.neighborChs
		%localKappaData.localKappaMatrix
		%deltaThetaAlphaMeanKappa
		%fds
		chanNum=bestLocalChanNum;
	end


        %localLFPdata=matfile(sprintf('C:/Users/tibin/spikeDynamicsAnalysisTibin/session3Mat/concatChan%d/lfpDecConcat.mat',chanNum));

	if(strcmp(subj,'MG49'))
        	localLFPdata=matfile(sprintf('/nfs/turbo/lsa-ojahmed/tibin/spikeDynamicsAnalysisTibin/session3Mat/concatChan%d/lfpDecConcat.mat',chanNum));
		%decFs=cell.Fs/15; %make sure this dec factor holds!!!!
		decFs=30000/15; %make sure this dec factor holds!!!!
		decLFP=localLFPdata.concatenatedLFP;
		return
	end

	%if(strcmp(subj,'thetaRat24_RSC_V1_HPC'))
	if(strcmp(subj,'FangChiData'))
		sourceData='/nfs/turbo/lsa-ojahmed/tibin/processedHumanData/FangChiData/sessionID-1';
		lfpData=load(fullfile(sourceData,'2018-02-13_14-00-20_Rat0024_Sleep_2HR_decLFP.mat'));
		decLFP=lfpData.decLFP;
		decFs=2000;
		return
	end

	%generic case
	chStr=getChStr(chanNum);
	lfpDataFilePath=getRegexFilePath(sprintf('/nfs/turbo/lsa-ojahmed/tibin/processedHumanData/%s/sessionID-%d/decimatedLFP-MatFiles/Fs2000Hz/',subj,sessionID),sprintf('DecLFP*_ch%s.mat',chStr));
	localLFPdata=matfile(lfpDataFilePath);
	decFs=localLFPdata.decFs;
	decLFP=localLFPdata.decLFP;
	

	%ns5FilesDir=sprintf('/nfs/turbo/lsa-ojahmed/tibin/processedHumanData-concatenated/MG63/05-103012-062_064_to_072',sub/j);

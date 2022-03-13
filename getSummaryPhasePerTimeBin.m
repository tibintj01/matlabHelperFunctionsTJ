function [summaryPhasePerCycle] = getSummaryPhasePerTimeBin(currFieldTimeCenters,currFieldAllTimes,currFieldSpikePhases,gaussKernelWidth)
%gets most common phase at each time bin center supplied for continuous
%times and phases

    useCircMean=1;
  useCircMean=0;

    phaseKernelDistrBinCenters=1:360;
    %gaussKernelWidth=30;
    
    if(~exist('gaussKernelWidth','var'))
     gaussKernelWidth=45;
     %%gaussKernelWidth=15;
    end

    numTimeCenters=length(currFieldTimeCenters);
    summaryPhasePerCycle=NaN(size(currFieldTimeCenters));
    %minNumPhasesForDistrEst=3;
       minNumPhasesForDistrEst=1;

    for ci=1:numTimeCenters
        currTimeCenter=currFieldTimeCenters(ci);

        currTimeBinIdxes=abs(currFieldAllTimes-currTimeCenter)<0.001;

        if(sum(currTimeBinIdxes)<minNumPhasesForDistrEst)
            continue
        end
        currTimeBinPhases=currFieldSpikePhases(currTimeBinIdxes);


        if(useCircMean)
            if(~isempty(currTimeBinPhases))
             summaryPhasePerCycle(ci)=circMeanDeg(currTimeBinPhases);
            end
        else
            vfObservations=currTimeBinPhases(:)*2*pi/360;
             vfPDFSamples=phaseKernelDistrBinCenters(:)*2*pi/360;
             vfDomain=[0 2*pi];
             fSigma=gaussKernelWidth*2*pi/360;
             [vfEstimate] = circ_ksdensity(vfObservations, vfPDFSamples, vfDomain, fSigma);

             currPdist=vfEstimate/sum(vfEstimate(:));
             %plot(phaseKernelDistrBinCenters,currPdist,'k-','LineWidth',3)

            [~,maxID]=max(currPdist);

            currBinSummaryPhase=phaseKernelDistrBinCenters(maxID);

            summaryPhasePerCycle(ci)=currBinSummaryPhase;
        end
    end
    
   
    


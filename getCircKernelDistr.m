function [currPdist,peakPhase] = getCircKernelDistr(phases,phaseEdges,gaussKernelWidth)
%

      if(~exist('phaseEdges','var'))
            phaseEdges=linspace(0,360,361);
        end
        
        if(~exist('gaussKernelWidth','var'))
            gaussKernelWidth=15;
        end
        phaseBinCenters=edgesToBins(phaseEdges);

            vfObservations=phases(:)*2*pi/360;
                 vfPDFSamples=phaseBinCenters(:)*2*pi/360;
                 vfDomain=[0 2*pi];
                 fSigma=gaussKernelWidth*2*pi/360;
                 [vfEstimate] = circ_ksdensity(vfObservations, vfPDFSamples, vfDomain, fSigma);
                 
                 currPdist=vfEstimate/sum(vfEstimate(:));
                 
                 
                 [maxVal,maxID]=max(currPdist);
            
            maxIDs=find(abs(currPdist-maxVal)<0.0001);
            
            
            peakPhase=nanmean(phaseBinCenters(maxIDs));

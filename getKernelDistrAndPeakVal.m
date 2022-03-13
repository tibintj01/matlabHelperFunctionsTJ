function [currPdist,peakVal] = getKernelDistrAndPeakVal(vals,gaussKernelSD,distrEdges,fH)
%

        distrBinCenters=edgesToBins(distrEdges);
        
        [kEst,xi] = ksdensity(vals,distrBinCenters,'Bandwidth',gaussKernelSD);
                %[kEst,xi] = ksdensity(vals,distrBinCenters);
                


         currPdist=kEst/sum(kEst(:));


         [maxVal,maxID]=max(currPdist);
            
            maxIDs=find(abs(currPdist-maxVal)<0.0001);
            
            
            peakVal=nanmean(distrBinCenters(maxIDs));
            
            if(exist('fH','var'))
                hold on
                plot(distrBinCenters,currPdist,'k-','LineWidth',3)
                hold on
                vline(peakVal,'k--',2)
                
            end

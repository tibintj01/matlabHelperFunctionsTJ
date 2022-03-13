function [xBinCenters,peakPhaseValues] = getBinnedCircPeakVals(x,phases,xEdges)
%computes circ average value of y in x edges specified 

gaussKernelWidthDeg=15;
    validIdxes=~isnan(x(:)) & ~isnan(phases(:));
    
    if(min(phases)>=0)
    distrPhaseEdges=linspace(0,360,360+1);
    else
         distrPhaseEdges=linspace(-180,180,360+1);
    end
    
    x=x(validIdxes);
    phases=phases(validIdxes);
    
    xBinCenters=edgesToBins(xEdges);
    xBinWidth=median(diff(xBinCenters));
    
    numXbins=length(xBinCenters);
    
    peakPhaseValues=NaN(numXbins,1);
    upperLimPerBin=NaN(numXbins,1);
    lowerLimPerBin=NaN(numXbins,1);
    for xi=1:numXbins
        
        currXbinStart=xBinCenters(xi)-xBinWidth/2;
         currXbinEnd=xBinCenters(xi)+xBinWidth/2;
         
         currBinIdxes=x>=currXbinStart & x<currXbinEnd;
         
         if(sum(currBinIdxes)==0)
             continue
         end
         
         if(xi==numXbins)
             currBinIdxes=x>=currXbinStart & x<=currXbinEnd;%include far edge in last bin
         end
         
         currBinPhases=phases(currBinIdxes);
         if(~isempty(currBinPhases))

            [~,currPeakPhase] = getCircKernelDistr(currBinPhases,distrPhaseEdges,gaussKernelWidthDeg);
                 
             peakPhaseValues(xi)=currPeakPhase;
         end
         

    end
    
    validIdxes=~isnan(peakPhaseValues(:)) & ~isnan(xBinCenters(:));

    
    
    peakPhaseValues=peakPhaseValues(validIdxes);
    xBinCenters=xBinCenters(validIdxes);

    
    xBinCenters=xBinCenters(:);
    peakPhaseValues=peakPhaseValues(:);

    

    
    
    

    


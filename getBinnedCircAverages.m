function [xBinCenters,phaseAvgValues,lowerLimPerBin,upperLimPerBin] = getBinnedCircAverages(x,phases,xEdges)
%computes circ average value of y in x edges specified 

    validIdxes=~isnan(x(:)) & ~isnan(phases(:));
    
    x=x(validIdxes);
    phases=phases(validIdxes);
    
    doNotModBy360=0;
    if(min(phases)<0)
        doNotModBy360=1;
    end
    xBinCenters=edgesToBins(xEdges);
    xBinWidth=median(diff(xBinCenters));
    
    numXbins=length(xBinCenters);
    
    phaseAvgValues=NaN(numXbins,1);
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
             if(doNotModBy360)
                 phaseAvgValues(xi)=circMeanDegNoMod360(currBinPhases);
             else
                phaseAvgValues(xi)=circMeanDeg(currBinPhases);
             end
         end
         
         [h mu ul ll]=circ_mtest(ang2rad(currBinPhases),ang2rad(phaseAvgValues(xi)));
         
         upperLim=rad2ang(ul);
         lowerLim=rad2ang(ll);
         
         if(~doNotModBy360 & isreal(upperLim) & isreal(lowerLim))
             upperLim=mod(upperLim,360);
             lowerLim=mod(lowerLim,360);
         end
         
        
         upperLimPerBin(xi)=upperLim;
         lowerLimPerBin(xi)=lowerLim;
        
    end
    
    validIdxes=~isnan(phaseAvgValues(:)) & ~isnan(xBinCenters(:));
    validIdxes=validIdxes & isreal(lowerLimPerBin) & isreal(lowerLimPerBin);
    
    
    phaseAvgValues=phaseAvgValues(validIdxes);
    xBinCenters=xBinCenters(validIdxes);
    upperLimPerBin=upperLimPerBin(validIdxes);
    
    lowerLimPerBin=lowerLimPerBin(validIdxes);
    
    xBinCenters=xBinCenters(:);
    phaseAvgValues=phaseAvgValues(:);
    lowerLimPerBin=lowerLimPerBin(:);
    upperLimPerBin=upperLimPerBin(:);
    

    
    
    

    


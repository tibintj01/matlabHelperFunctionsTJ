function [xBinCenters,yAvgValues, ySEMPerBin,yStdDevPerBin ] = getBinnedAverages(x,y,edges)
%computes average value of y in x edges specified 


    minNumPts=0;
    
    validOriginalIdxes=~isnan(x(:)) & ~isnan(y(:)) & x>=min(edges) & x<=max(edges);
    
    x=x(validOriginalIdxes);
    y=y(validOriginalIdxes);
    
    if(isempty(x) || isempty(y))
        xBinCenters=NaN;
        yAvgValues=NaN;
        return
    end
    
    
    [~,~,loc]=histcounts(x,edges);
    
    usedXBinIdxes=unique(loc);
    
    numValsPerBin=accumarray(loc(:),1);
    yAvgValues = accumarray(loc(:),y(:))./numValsPerBin;
    
    ySEMPerBin=NaN(size(yAvgValues));
    yStdDevPerBin=NaN(size(yAvgValues));
    
    for xi=1:max(loc)
        currXBinIdxes=(loc==xi);
        currXbinYvals=y(currXBinIdxes);
        
        ySEMPerBin(xi)=getSEMacrossRows(currXbinYvals(:));
        
        yStdDevPerBin(xi)=nanstd(currXbinYvals(:));
    end
    
    xBinCenters = 0.5*(edges(1:end-1)+edges(2:end));
    
    xBinCenters=xBinCenters(usedXBinIdxes);
    yAvgValues=yAvgValues(usedXBinIdxes);
    
    lowCountIdxes=numValsPerBin<minNumPts;
    
    xBinCenters(lowCountIdxes)=NaN;
    yAvgValues(lowCountIdxes)=NaN;
    
    
    validAvgIdxes=~isnan(xBinCenters(:)) & ~isnan(yAvgValues(:));
    
    xBinCenters=xBinCenters(validAvgIdxes);
    yAvgValues=yAvgValues(validAvgIdxes);
    
    ySEMPerBin=ySEMPerBin(validAvgIdxes);
    yStdDevPerBin=yStdDevPerBin(validAvgIdxes);
    


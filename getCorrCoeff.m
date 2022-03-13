function R=getCorrCoeff(x,y)
%UNTITLED Summary of this function goes here

    x=x(:);
    y=y(:);
    goodIdxes=~isnan(x) & ~isnan(y);
    
    rMat=corrcoef(x(goodIdxes),y(goodIdxes));
    
    if(length(rMat)==1)
        R=NaN;
    else
        
        R=rMat(1,2);
    end


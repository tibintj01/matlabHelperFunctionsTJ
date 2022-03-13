function [zp] = getLinearCrossSection(z,startRowCol,endRowCol,numPts)
%returns z along parameterized line specified by start and end coordinates
    %p=linspace(0,1,round(numPts));
    if(~exist('numPts','var'))
        numPts=100;
    end
    
    rowSeq=linspace(startRowCol(1), endRowCol(1),numPts);
    colSeq=linspace(startRowCol(2), endRowCol(2),numPts);
    
    for pi=1:numPts
         currRow=round(rowSeq(pi));
         currCol=round(colSeq(pi));
         zp(pi)=z(currRow,currCol);
    end
end


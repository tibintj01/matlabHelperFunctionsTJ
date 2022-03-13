function [surrogateRelEntropies]=testBootstrapKLDivFull(pVals,qVals,edges,Nbootstrap)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
disp('getting surrogate relative entropies....')
tic

isCirc=0;
showPlots=1;

    nPvals=length(pVals);
    nQvals=length(qVals);
    
    combinedPopVals=[pVals(:); qVals(:)];
    
    surrogateRelEntropies=NaN(Nbootstrap,1);
    for i=1:Nbootstrap
        if(mod(i,1000)==0)
            disp(sprintf('finished %d surrogate values...',i))
        end
        pSurrogateVals=randsample(combinedPopVals,nPvals);
        qSurrogateVals=randsample(combinedPopVals,nQvals);
        
        [pSurrogate] = getProbDist(pSurrogateVals,edges,showPlots,isCirc);
        [qSurrogate] = getProbDist(qSurrogateVals,edges,showPlots,isCirc);
        
        surrogateRelEntropies(i)=getRelativeEntropy(pSurrogate,qSurrogate);
    end
    
    toc
    
   


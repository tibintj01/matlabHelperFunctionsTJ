function [surrogateMeanDiffs,originalMeanDiff]=getBootstrapTwoGpMeanCircDiff(gp1Vals,gp2Vals,Nbootstrap)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
disp('getting surrogate angular differences....')
tic

isCirc=0;
showPlots=1;

    n1vals=length(gp1Vals);
    n2vals=length(gp2Vals);
    
    gp2MinusGp1=NaN(n1vals,1);
    
    for i=1:n1vals
        gp2MinusGp1(i)=angdiffDeg([gp1Vals(i),gp2Vals(i)]);
        
    end
    
    originalMeanDiff=nanmean(gp2MinusGp1);
    
    combinedPopVals=[gp1Vals(:); gp2Vals(:)];
    
    surrogateMeanDiffs=NaN(Nbootstrap,1);
    for i=1:Nbootstrap
        if(mod(i,1000)==0)
            disp(sprintf('finished %d surrogate values...',i))
        end
        gp1SurrogateVals=randsample(combinedPopVals,n1vals);
        gp2SurrogateVals=randsample(combinedPopVals,n2vals);
        
        gp2MinusGp1=NaN(n1vals,1);
    
        for j=1:n1vals
            gp2MinusGp1(j)=angdiffDeg([gp1SurrogateVals(j),gp2SurrogateVals(j)]);

        end
        
        currSurrogateMeanDiff=nanmean(gp2MinusGp1);

        surrogateMeanDiffs(i)=currSurrogateMeanDiff;
    end
    
    toc
    
   


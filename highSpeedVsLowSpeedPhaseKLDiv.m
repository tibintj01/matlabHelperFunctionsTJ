clear all;  clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%statistically test self-similarity of distributions up to slide
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data=load('lowSpeedVsHighSpeedPhaseProbDist.mat')

testTimeInFieldSlideSimilarity=1;
testTimeInFieldSlideSimilarity=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(testTimeInFieldSlideSimilarity)
    allFieldMeanValsHighSpeed=data.masterMeanTimeToHighSpeeds;
    allFieldMeanValsLowSpeed=data.masterMeanTimeToLowSpeeds;
else
    allFieldMeanValsHighSpeed=data.allFieldMeanPhasesHighSpeed;
    allFieldMeanValsLowSpeed=data.allFieldMeanPhasesLowSpeed;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%recalculate distributions using own edges
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%phaseEdgesDeg=data.phaseEdgesDeg;

if(testTimeInFieldSlideSimilarity)
  
    numBins=200;
      valEdges=linspace(0,5,numBins+1);
      isCirc=0;
else
      numBins=50;
    numBins=40;
    %numBins=30;
    valEdges=linspace(0,360,numBins+1);
    isCirc=1;
      
end

pSlow=getProbDist(allFieldMeanValsLowSpeed,valEdges,1,isCirc);
pFast=getProbDist(allFieldMeanValsHighSpeed,valEdges,1,isCirc);
%{
pSlow=data.pPhaseSlow;
pFast=data.pPhaseFast;
%}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get distributions zero mean data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(isCirc)
    meanValHigh=circMeanDeg(allFieldMeanValsHighSpeed);
    meanValLow=circMeanDeg(allFieldMeanValsLowSpeed);
else
    meanValHigh=nanmean(allFieldMeanValsHighSpeed);
    meanValLow=nanmean(allFieldMeanValsLowSpeed);
end


zeroedValsHighSpeed=mod(allFieldMeanValsHighSpeed-meanValHigh,360);
zeroedValsLowSpeed=mod(allFieldMeanValsLowSpeed-meanValLow,360);


pSlowZeroed=getProbDist(zeroedValsLowSpeed,valEdges,1,isCirc);
pFastZeroed=getProbDist(zeroedValsHighSpeed,valEdges,1,isCirc);

autoArrangeFigures


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get relative entropy of original data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[DslowFast] = getRelativeEntropy(pSlow,pFast)
[DslowFastZeroed] = getRelativeEntropy(pSlowZeroed,pFastZeroed)

N_bootstrap=100000;
%Nbootstrap=50000;
Nbootstrap=10000;

[surrogateRelEntropiesOriginal]=getBootstrappedRelativeEntropies(allFieldMeanValsLowSpeed,allFieldMeanValsHighSpeed,valEdges,N_bootstrap);

%%
figure
%histogram(surrogateRelEntropiesOriginal,100);
[Nbootstrap,bootstrapEdges]=histcounts(surrogateRelEntropiesOriginal,100);
bootstrapBins=edgesToBins(bootstrapEdges);

plot(bootstrapBins,Nbootstrap,'k-','LineWidth',2)


hold on



pHoriginal=plot([DslowFast DslowFast],ylim,'k-','LineWidth',3)
pHzeroed=plot([DslowFastZeroed DslowFastZeroed],ylim,'k--','LineWidth',3)

[pValOriginal]=getBootstrapPVal(Nbootstrap,bootstrapBins,DslowFast);
[pValZeroed]=getBootstrapPVal(Nbootstrap,bootstrapBins,DslowFastZeroed);

title(sprintf('Bootstrapped KL divergence, N_{bootstrap}=%d, Original data, p=%.4f; Mean-zeroed data, p=%.4f',N_bootstrap,pValOriginal,pValZeroed))

legend([pHoriginal, pHzeroed], {'original data', 'mean-zeroed data'})
setFigFontTo(32)
maxFigHalfWidth
box off
legend boxoff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get relative entropy of zero mean distributions - same as nonzero mean,
%entropy is not sensitive to mean shift!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{

[surrogateRelEntropiesZeroed]=getBootstrappedRelativeEntropies(zeroedPhasesLowSpeed,zeroedPhasesHighSpeed,phaseEdgesDeg,Nbootstrap);

figure
histogram(surrogateRelEntropiesZeroed,100);
%}




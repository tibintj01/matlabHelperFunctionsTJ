function [phaseBinCentersPadded,currPdistPadded,padLength]=getCircPaddedDistr(phaseBinCenters,currPdist,upSampleFact)

if(~exist('upSampleFact','var'))
    upSampleFact=1;
end



midIdxBottom=ceil(length(phaseBinCenters)/2);
midIdxTop=floor(length(phaseBinCenters)/2);

bottomPhases=phaseBinCenters(midIdxBottom:end) - 360;
topPhases=phaseBinCenters(1:midIdxTop) +360;

bottomProbs=currPdist(midIdxBottom:end);
topProbs=currPdist(1:midIdxTop);

phaseBinCentersPadded=[bottomPhases(:); phaseBinCenters(:); topPhases(:)];

currPdistPadded=[bottomProbs(:) ;currPdist(:); topProbs(:)];

padLength=length(bottomPhases)*upSampleFact;


numBinsOriginal=length(phaseBinCentersPadded);
numBinsUpsampled=numBinsOriginal*upSampleFact;

upsampledBins=linspace(1,numBinsOriginal,numBinsUpsampled);

phaseBinCentersUpSampled=interp1(1:numBinsOriginal,phaseBinCentersPadded,upsampledBins);
currPdistUpSampled=interp1(1:numBinsOriginal,currPdistPadded,upsampledBins);

phaseBinCentersPadded=phaseBinCentersUpSampled;
currPdistPadded=currPdistUpSampled;

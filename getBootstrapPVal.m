function [pVal]=getBootstrapPVal(Nbootstrap,bootstrapBins,originalVal)
%get p value given surrogate distribution and original value
%(area under distirubiton to the right of original value)

aboveIdxes=bootstrapBins>=originalVal;

pBootstrap=Nbootstrap/sum(Nbootstrap(:));

pVal=sum(pBootstrap(aboveIdxes));

if(pVal>0.5)
    pVal=1-pVal;
end


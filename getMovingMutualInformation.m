function [localPosPhaseMIindex,posBinCenters]=getMovingMutualInformation(spikePositions, spikePhases, trackLengthM )
%given spike positions and phases, checks each position bin for significant
%phase modulation in the form of entropy relative to max entropy, i.e. deviation from uniform

useJointToEstEntropy=1;


%posBinWidth=0.02; %m
%posBinWidth=0.03; %m
%posBinWidth=0.05; %m
posBinWidth=0.025; %m
%posBinWidth=0.10; %m
posBinEdges=0:posBinWidth:trackLengthM;
posBinCenters=edgesToBins(posBinEdges);

numPosBins=length(posBinCenters);

phaseBinWidth=15;
phaseBinWidth=20;
%phaseBinWidth=30;
%phaseBinWidth=10;
%phaseBinWidth=5;
phaseEdges=0:phaseBinWidth:360;
phaseBinCenters=edgesToBins(phaseEdges);

numPhaseBins=length(phaseBinCenters);

posWindSize=0.1; %m
posWindSize=0.2; %m
posWindSize=0.4; %m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%generated 2d joint distribution of phase and position
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



isCirc=1;
showPlots=0;

%minNumSpikesInBin=3;
minNumSpikesInBin=10;
minNumSpikesInBin=5;
minNumSpikesInBin=3;
minNumSpikesInBin=20;

localPosPhaseMI=NaN(size(posBinCenters));
localPosPhaseMIindex=NaN(size(posBinCenters));
%phaseEntropyIdxPerPos=NaN(size(posBinCenters));

%minNumSpikesPerPhaseDistr=10;
%{
jointProbThresh=0.1;
jointProbThresh=0.01;
jointProbThresh=0.0075;
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%loop through position bins and get phase distribution per bin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for pi=1:length(posBinCenters)
    
    currPosStart=posBinCenters(pi)-posWindSize/2;
    currPosEnd=posBinCenters(pi)+posWindSize/2;
    
    currPosBinSpikeIdxes=spikePositions>currPosStart & spikePositions<currPosEnd;
    
    if(sum(currPosBinSpikeIdxes)<minNumSpikesInBin)
            continue
    end
    
    currBinSpikesPhases=spikePhases(currPosBinSpikeIdxes);
    currBinSpikesPositions=spikePositions(currPosBinSpikeIdxes);
    posBinEdgesLocal=currPosStart:posBinWidth:currPosEnd;
    
    posBinCentersLocal=edgesToBins(posBinEdgesLocal);
    
    if(showPlots)
        fH=figure(19)

        subplot(2,1,1)
        pPosPhase=getJointDistr(currBinSpikesPositions,currBinSpikesPhases,posBinEdgesLocal,phaseEdges,fH);

    else
       pPosPhase=getJointDistr(currBinSpikesPositions,currBinSpikesPhases,posBinEdgesLocal,phaseEdges);

    end
     
    pPos=getProbDist(currBinSpikesPositions,posBinEdgesLocal,0,0);
    pPhase=getProbDist(currBinSpikesPhases,phaseEdges,0,1);
    
    independenceDistr=NaN(length(pPos),length(pPhase));
    for ii=1:length(pPos)
        for jj=1:length(pPhase)
            independenceDistr(ii,jj)=pPos(ii)*pPhase(jj);
        end
    end
    
    if(showPlots)
        subplot(2,1,2)
        omarPcolor(posBinCentersLocal,phaseBinCenters,independenceDistr',fH)
        cb=colorbar
        colormap(gca,jet)
        ylabel(cb,'probability if independent')
    end
    localPosPhaseMI(pi)=getRelativeEntropy2D(pPosPhase,independenceDistr);
    maxPossibleMI=log2(length(pPosPhase(:)));

    localPosPhaseMIindex(pi)=(localPosPhaseMI(pi))/maxPossibleMI;

    %close(fH)
end


%close all



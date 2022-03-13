function [phaseEntropyIdxPerPos,posBinCenters,sigThresh]=getPhaseDistEntropyByPos(spikePositions, spikePhases, trackLengthM )
%given spike positions and phases, checks each position bin for significant
%phase modulation in the form of entropy relative to max entropy, i.e. deviation from uniform



useJointToEstEntropy=1;
%useJointToEstEntropy=0;
%Nshuffles=100;

posBinWidth=0.02; %m
posBinWidth=0.03; %m
posBinWidth=0.05; %m
%posBinWidth=0.025; %m
%posBinWidth=0.10; %m
posBinEdges=0:posBinWidth:trackLengthM;
posBinCenters=edgesToBins(posBinEdges);

numPosBins=length(posBinCenters);

phaseBinWidth=15;
%phaseBinWidth=20;
%phaseBinWidth=30;
%phaseBinWidth=10;
%phaseBinWidth=5;
phaseEdges=0:phaseBinWidth:360;
phaseBinCenters=edgesToBins(phaseEdges);

numPhaseBins=length(phaseBinCenters);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%generated 2d joint distribution of phase and position
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(useJointToEstEntropy)
    pPosPhase=getJointDistr(spikePositions,spikePhases,posBinEdges,phaseEdges);
    %pPos=getProbDist(spikePositions,posBinEdges,0,0);
    %pPhase=getProbDist(spikePhases,phaseEdges,0,1);

end
maxPossibleEntropy=log2(numPhaseBins);

isCirc=1;
showPlots=0;

%minNumSpikesInBin=3;
minNumSpikesInBin=10;
minNumSpikesInBin=5;
minNumSpikesInBin=3;

phaseHPerPos=NaN(size(posBinCenters));
phaseEntropyIdxPerPos=NaN(size(posBinCenters));

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
    
    currPosStart=posBinCenters(pi)-posBinWidth/2;
    currPosEnd=posBinCenters(pi)+posBinWidth/2;
    
    currPosBinSpikeIdxes=spikePositions>currPosStart & spikePositions<currPosEnd;
    
    currPosBinSpikesPhases=spikePhases(currPosBinSpikeIdxes);
    
      if(length(currPosBinSpikesPhases)<minNumSpikesInBin)
            continue
        end

    if(~useJointToEstEntropy)
        
            [currPosPhaseDist] = getProbDist(currPosBinSpikesPhases,phaseEdges,showPlots,isCirc);
    else
            
            totalProbThisPos=sum(smooth(pPosPhase(pi,:)));
            currPosPhaseDist=smooth(pPosPhase(pi,:))/totalProbThisPos;
            
            %{
            figure(111)
            subplot(2,1,1)
            plot(spikePositions, spikePhases,'k.')
            ylim([0 360])
            subplot(2,1,2)
            plot(phaseBinCenters,currPosPhaseDist)
            xlim([0 360])
           hold on
           title(sprintf('Pos: %.2f, H Idx: %.2f',posBinCenters(pi),(maxPossibleEntropy-getEntropy(currPosPhaseDist))/maxPossibleEntropy))
            %}
    
    end
            phaseHPerPos(pi)=getEntropy(currPosPhaseDist);
            phaseEntropyIdxPerPos(pi)=(maxPossibleEntropy-phaseHPerPos(pi))/maxPossibleEntropy; %map to 0 to 1 scale of phase-space coding

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %try uniform distribution at each position to get significance
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
    
end
close all



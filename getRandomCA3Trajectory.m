function [ca3Trajectory] = getRandomCA3Trajectory(nCA3,rndSeed,seqBased)
%given number of CA3 cells and random initialization seed (EC3)
%returns CA3 trajectory generated as follows
%CA3 sequence simulated as random selection of 7 cells
%updated by increasing sequence by one (phase precession)
%and replacing maxed-out cell by one random CA3 cell
%and repeating ad infinitum; 
%-->CA1 recruitment stats? Field number distribution?
numIterations=5000;
numIterations=7500;
numIterations=1500;
numIterations=2000;
numIterations=7000;
numIterations=6000;
numIterations=2000;
%numIterations=1000;
numIterations=1500;
rng(rndSeed)

initial7=randsample(nCA3,7,0);


updatedCA3seq=initial7;

ca3Trajectory=NaN(numIterations,nCA3);
%if(seqBased)
    for i=1:numIterations
        %ca3Trajectory(i,:)=updatedCA3seq;
        ca3Trajectory(i,:)=listToActivationVec(updatedCA3seq,nCA3);
        updatedCA3seq=circshift(updatedCA3seq,1);
        %find unit not in current line up to add
        newUnit=randsample(nCA3,1);
        while(min(abs(updatedCA3seq-newUnit))==0)
            newUnit=randsample(nCA3,1);
        end
        updatedCA3seq(1)=newUnit;
    end
    
    %if(~seqBased)
    %    ca3Trajectory=ca3Trajectory(randperm(numIterations),randperm(nCA3));
    %end
%else
%    for i=1:numIterations
%        %ca3Trajectory(i,:)=listToActivationVec(repelem(randsample(nCA3,1),7),nCA3);
%        ca3Trajectory(i,:)=listToActivationVec(randsample(nCA3,7,0),nCA3);
%    end
%end





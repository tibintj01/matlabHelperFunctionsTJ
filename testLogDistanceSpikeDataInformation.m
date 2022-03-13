%%
%close all
%data=load('allFieldSpaceTimePhaseTriplets.mat');
data=load('allFieldSpaceTimePhaseTripletsAndExcessJoint.mat');
allFieldSpaceTimePhaseTriplets=data.allFieldSpaceTimePhaseTriplets;

spaceVals=allFieldSpaceTimePhaseTriplets(:,1);

spikePhases=allFieldSpaceTimePhaseTriplets(:,3);

%figure;

%plot(spaceVals,spikePhases,'k.')

%%

%numSpaceBins=30;
%numPhaseBins=30;

numSpaceBins=max(spaceVals)+1; % for discretize (0 ->1)
%numSpaceBins=16;
%numPhaseBins=90;

%numPhaseBins=20;
%numPhaseBins=30;
%numSpaceBins=15;
%numPhaseBins=10;

permIndices=randperm(length(spaceVals));
permIndices=1:length(spaceVals);

%IspacePhase=getDiscreteMutualInfo(log(spaceVals),spikePhases,numSpaceBins,numPhaseBins);
numPhaseBins=45;
numPhaseBins=30;

phaseBinNums=3:2:45;
%phaseBinNums=45:45;
phaseBinNums=20:20;
phaseBinNums=30:30;
%phaseBinNums=3:3:90;
%phaseBinNums=45:45;
%phaseBinNums=16:16;
%phaseBinNums=3:2:180;
numPhaseBinSizes=length(phaseBinNums);

ItoPhaseEntropy=zeros(numPhaseBinSizes,1);
IspacePhaseAll=zeros(numPhaseBinSizes,1);
HphaseAll=zeros(numPhaseBinSizes,1);
normSpaceVals=(spaceVals(permIndices))/max(spaceVals);

scaleFact=10e5;
scaleFact=1;
rescaleSpaceVals=(normSpaceVals)*scaleFact;
%logSpaceVals=log(scaleFact-rescaleSpaceVals);

disp('here')
numShuffles=1000;
%numShuffles=1;
IspacePhaseShuffles=NaN(numShuffles,1);

for pbi=1:numPhaseBinSizes
    numPhaseBins=phaseBinNums(pbi);
    numSpaceBins=phaseBinNums(pbi);
    for shufflei=1:numShuffles
        if(shufflei==1)
            permIndices=1:length(spaceVals);
	   showPlotShuffle=1;
        else
            permIndices=randperm(length(spaceVals));
	   showPlotShuffle=0;
        end
        %[IspacePhase,Hspace,Hphase]=getDiscreteMutualInfo(normSpaceVals,(spikePhases),numSpaceBins,numPhaseBins,'Dist in field','Theta phase (deg)')
       %[IspacePhase,Hspace,Hphase]=getDiscreteMutualInfo((rescaleSpaceVals),(spikePhases),numSpaceBins,numPhaseBins,'Dist in field','Theta phase (deg)')
        [IspacePhase,Hspace,Hphase]=getDiscreteMutualInfo((rescaleSpaceVals(permIndices)),(spikePhases),numSpaceBins,numPhaseBins,'Dist in field','Theta phase (deg)',showPlotShuffle);
        
        IspacePhaseShuffles(shufflei)=IspacePhase;
    end
        ItoPhaseEntropy(pbi)=IspacePhase/Hphase;
        IspacePhaseAll(pbi)=IspacePhase;
        HphaseAll(pbi)=Hphase; 
end

figure
histogram(IspacePhaseShuffles)
hold on

plot([IspacePhaseShuffles(1) IspacePhaseShuffles(1)],ylim,'k-','LineWidth',5)


%{

figure;

plot(360./phaseBinNums,ItoPhaseEntropy)
xlabel('Theta phase bin size (deg)')
ylabel('Space-phase mutual information to phase entropy ratio')

figure
plot(360./phaseBinNums,HphaseAll)
xlabel('Theta phase bin size (deg)')
ylabel('Phase entropy')

figure;
plot(360./phaseBinNums,IspacePhaseAll)
xlabel('Theta phase bin size (deg)')
ylabel('Space-phase mutual information')

figure;
plot(phaseBinNums,IspacePhaseAll)
xlabel('number of bins (phase and space)')
ylabel('Space-phase mutual information')
%}


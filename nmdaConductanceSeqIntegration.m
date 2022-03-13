

dt=0.01; %msec
simTime=200 %msec
numSteps=round(simTime/dt);

tAxis=dt:dt:simTime;

numCompartments=7;
tau2=10; %msec
tau1s=linspace(20,50,numCompartments);

tau1s=fliplr(tau1s);

spatiotemporalSpikePattern=zeros(numCompartments,numSteps);

gapTimeIdx=round(20/dt);

for i=1:numCompartments
	spikeIdx=gapTimeIdx*(i-1)+1;
	spatiotemporalSpikePattern(i,spikeIdx)=1;
end

gsynMax=1;
gsyn=zeros(numCompartments,numSteps);
weights=gsynMax*ones(numCompartments,1);

for step=1:numSteps-1
	%loop across dendritic compartments (towards soma)
	%for compNum=1:numCompartments
	currentSpatialSpikePattern=spatiotemporalSpikePattern(:,step);

	currentCompsWithSpikes=find(currentSpatialSpikePattern);
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%for each spike-compartment pair, adds kernel to conductance variable
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	for i=1:length(currentCompsWithSpikes)
		compNum=currentCompsWithSpikes(i);
		temporalCoeff=tau2*tau1s(compNum)/(tau1s(compNum)-tau2);

		currSpikeTime=step*dt;
		weight=weights(compNum);

		synEndStep=step+1+round(300/dt);
		if(synEndStep>numSteps)
			synEndStep=numSteps;
		end
		synCurrentIdxes=(step+1):synEndStep;

		gsyn(compNum,synCurrentIdxes)=weight*temporalCoeff*exp(-(dt*synCurrentIdxes-currSpikeTime)/tau1s(compNum))-exp(-(dt*synCurrentIdxes-currSpikeTime)/tau2);
	
	end

end

figure
plot(tAxis,gsyn)

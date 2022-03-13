function [] =seqDrivenPlaceNetworkPar(speedFact,noiseConditionID)

assemblyN=50;
normDIall=NaN(assemblyN,1);
rawDIall=NaN(assemblyN,1);


	%clearvars -except assemblyN noiseConditionID speedFact normDIall rawDIall
	%close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%speedFact=0.25;
%speedFact=1;
%speedFact=0.5;
N=1000;
%N=200;
minNumVisitedCells=round(N*0.95);
%N=500;
%N=200;
%N=100;
%N=100;
%thresholdPrecessionStepRatio=1;
%thresholdPrecessionStepRatio=0.5;
thresholdPrecessionStepRatio=2; %makes recurrent place fields less likely?; no winner take all seems to dominate
%thresholdPrecessionStepRatio=2.1; %makes recurrent place fields less likely?; no winner take all seems to dominate
%thresholdPrecessionStepRatio=2.5; %makes recurrent place fields less likely?; no winner take all seems to dominate
%thresholdPrecessionStepRatio=3; %makes recurrent place fields less likely?; no winner take all seems to dominate
detectAssemblies=0;
plotDendriticDrive=0;
%plotDendriticDrive=1;

saveActiveSynapseMovie=0;

cellCoverageFig=figure
%TEST LOCAL COPY CHANGE

filePath=sprintf('../processedData/seqDrivenPlaceNetworkMeasurements%.4f_N_%d_NoiseCondition_%d.mat',speedFact,N,noiseConditionID);

disp(sprintf('computing for noise condition id %d.........',noiseConditionID))
tic
	rerun=1;
	%rerun=0;
	%if(rerun)
	%	if(exist(filePath,'file'))
	%		delete(filePath)	
	%	end
	%end

	if(rerun || ~exist(filePath,'file'))
		%close all
		disp('new')
		%rng(2)
		%rng(3)
		rngVal=1;
		rng(rngVal)

		%N=7000;
		%N=7000;
		%N=700;
		%N=5000;
		%N=1000;
		%N=100;
		%N=500;
		%speedFact=0.5;
		%speedFact=1;
		%speedFact=1.2;
		%speedFact=2;
		%N=1000;
		numCompartmentsPerNeuron=7;
		%assemblyN=30;
		%assemblyN=30;
		%assemblyN=100;
		%assemblyN=300;
		%assemblyN=400;
		%assemblyN=round(N/3);
		%assemblyN=500;
		%localConnectivityRadius=500;
		localConnectivityRadius=assemblyN/2-1;
		%assemblyN=50;
		%assemblyN=25;
		%startOffsetNeuron=400;
		%startOffsetNeuron=0;
		%startOffsetNeuron=round(N/3);
		%startOffsetNeuron=20;
		startOffsetNeuron=0;

		%usePointConnectivity=1;
		usePointConnectivity=0;

		useLocalConnectivity=0;
		noSelfConnections=1;

		sequenceMatchDriven=1;
		%sequenceMatchDriven=0;


		%numNoiseConditions=30;
		numNoiseConditions=assemblyN;
		frozenUniformNoiseConditions=rand([numNoiseConditions assemblyN]);

		%numDendrites=5040;
		%numDendrites=7000;
		%7 neural subpopulation units
		%N=10;
		%N=5040;
		%N=20;
		%N=50;
		%N=100;
		%N=700;
		%N=140;
		%numCycles=100;
		%numCycles=101;
		numCycles=1000;
		activitySeq=(1:7);
		%networkExcStrength=0.5;
		%networkExcStrength=0.101;
		if(usePointConnectivity)
			networkExcStrength=1;
			connectivityStr='point_connectivity';
		else
			%networkExcStrength=1.5;
			%networkExcStrength=2;
			%networkExcStrength=1;
			networkExcStrength=0.5;
			
			%networkExcStrength=0.3;
			%networkExcStrength=0.1;
			%networkExcStrength=0.075;
			%networkExcStrength=0.0625;
			%networkExcStrength=0.05;
			%networkExcStrength=0.05;
			connectivityStr='extended_dendrite_connectivity';
		end
		%networkExcStrength=0.05;
		activitySeq=activitySeq/norm(activitySeq);

		%maxCycleDisp=7;
		%maxCycleDisp=3;
		%maxCycleDisp=5;
		%maxCycleDisp=15;
		%maxCycleDisp=30;
		%maxCycleDisp=100;
		maxCycleDisp=2000;

		neuralActivityVectorMovie=NaN(N,numCycles);
		functionalConnectivityMovie=NaN(N,N,numCycles);
		
		dendriticActivationOverTime=NaN(N,numCycles);
		activityDistOverTime=NaN(N,numCycles);
		
		%useBackRamp=1;
		useBackRamp=0;

		%useFlatRamp=1;
		useFlatRamp=0;

		%cycleCMap=copper(maxCycleDisp)
		cycleCMap=parula(maxCycleDisp)

		%activityPrecessionPerCycle=max(diff(activitySeq));
		%maxActivity=max(activitySeq);
		%minActivity=min(activitySeq);

		%spikeThreshold=minActivity+(maxActivity-minActivity)/2;
		%spikeThreshold=minActivity+(maxActivity-minActivity)/4;

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%set initial neural activity state vector
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%a=zeros(N,numCycles);
		a=zeros(N,numCycles);

		%activitySeq=(1:7)/7;
		%startSingleNeurons=1;
		startSingleNeurons=0;
		if(startSingleNeurons)
			startingContextSequence=activitySeq;
			%startingNeuronIDs=randsample(N,7);
			%startingNeuronIDs=1:7;
			startingNeuronIDs=(N-6):N;
			a(startingNeuronIDs,1)=startingContextSequence;
		else
			smoothRamp=1;
			
			if(smoothRamp==1)
				%smooth ramp
				%startNeuron=N-7*assemblyN-startOffsetNeuron;
				%endNeuron=min(N,startNeuron+7*assemblyN)-startOffsetNeuron;
				startNeuron=N-assemblyN-startOffsetNeuron;
				%endNeuron=min(N,startNeuron+assemblyN)-startOffsetNeuron;
				endNeuron=min(N,startNeuron+assemblyN)-startOffsetNeuron-1;
				for n=startNeuron:endNeuron
					if(useFlatRamp)
						%a(n,1)=startNeuron+rand(1)*20;
						a(n,1)=endNeuron-0.6*n;
					elseif(useBackRamp)
						a(n,1)=endNeuron-n;
					else
						a(n,1)=n-startNeuron + assemblyN/10*frozenUniformNoiseConditions(noiseConditionID,n-startNeuron+1);
					end
					%a(n,1)=(startNeuron+endNeuron)/2-startNeuron;
				end
				%a(startNeuron:endNeuron,1)=a(startNeuron:endNeuron,1)/norm(a(startNeuron:endNeuron,1))*10;
				a(startNeuron:endNeuron,1)=a(startNeuron:endNeuron,1)/norm(a(startNeuron:endNeuron,1))*5;
			else
				%staircase ramp
				for i=1:7
					%startNeuron=(i-1)*floor(N/7)+1;
					%endNeuron=min(N,startNeuron+floor(N/7));
					startNeuron=N-(7-i+1)*assemblyN;
					endNeuron=min(N,startNeuron+assemblyN);
					a(startNeuron:endNeuron,1)=activitySeq(i); 
					%a(startNeuron:endNeuron,1)=activitySeq(8-i); %dies out without propagating
					%a(startNeuron:endNeuron,1)=randsample(activitySeq,1);
					%a(startNeuron:endNeuron,1)=2;
					%(startNeuron:endNeuron,1)=activitySeq(end-i+1);
				end
			end
		end
		maxActivity=max(a(:,1));
		minActivity=min(a(:,1));
		baselinePrecessionStepSize=(maxActivity-minActivity)/7;
		%activityPrecessionPerCycle=(maxActivity-minActivity)/7*speedFact;
		%activityPrecessionPerCycle=baselinePrecessionStepSize*speedFact;
		activityPrecessionPerCycle=baselinePrecessionStepSize*speedFact*0.625/2;

		if(activityPrecessionPerCycle==0)
			activityPrecessionPerCycle=0.1;
		end

		figure; plot(a(:,1))
		xlim([startNeuron endNeuron])
		xlabel('cell no.')
		ylabel('Activity level/Theta phase')
		title(sprintf('Initial pattern (noise condition %d)',noiseConditionID))
		setFigFontTo(16)	
		saveas(gcf,sprintf('../figures/initialPattern_%s_noiseID_%d.tif',connectivityStr,noiseConditionID))
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%reshuffle initial sequence to ensure orthogonality
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		reshuffledIDs=randperm(N);
		a=a(reshuffledIDs,:);

		%spikeThreshold=minActivity+(maxActivity-minActivity)/2;
		%spikeThreshold=minActivity+(maxActivity-minActivity)/4;
		%spikeThreshold=minActivity+(maxActivity-minActivity)/7;
		spikeThreshold=minActivity+baselinePrecessionStepSize*thresholdPrecessionStepRatio;
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%build dendritic sequence detector network
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%seqConnMatrix=NaN(N,N);
		%postSynN=10*N;
		%postSynN=6000; 
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%assuming functional network is highly divergent in this functional context
		%i.e. each neuron has more output neurons/synchronous populations (~50) than input neurons/populations (7)
		%or perhaps "types" of functional input neurons as defined by template order and gamma
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%postSynA=zeros(postSynN,numCycles);
		%seqConnMatrix=zeros(postSynN,N);
		%seqConnMatrix=zeros(N,7);
		%seqConnMatrix=zeros(numDendrites,N);
		%seqConnMatrix=zeros(N,numCompartmentsPerNeuron);
		seqConnMatrix=zeros(N,N);

		%postSynTemplateR=NaN(postSynN,1);
		%postSynTemplateR=NaN(N,1);

		allSequences=perms(1:7);
		%allSequences=perms(activitySeq);
		%for postSynID=1:postSynN
		%for postSynID=1:N
		for preSynID=1:N
		%for postSynID=1:numDendrites
			%enough divergence to make expected convergence to be 7 for this N
			%divNum=ceil(N/7);
			%divNum=ceil(N/12);
			divNum=ceil(N/15);
			%divNum=ceil(N/17);
			%divNum=ceil(N/18);
			%divNum=ceil(N/20);
			%divNum=ceil(N/25);
			%divNum=ceil(N/50);
			%divNum=ceil(N/100);
			%divNum=ceil(N/5);
			%divNum=ceil(N/4);
			%divNum=ceil(N);

			%divNum=ceil(N/3);
			%divNum=ceil(N/4);
			%divNum=ceil(N/5);
			%divNum=ceil(N/2.5);
			%divNum=ceil(N/2);
			%divNum=ceil(N/14);
			%divNum=ceil(N/14/2);
			%presynNeuronIDs=randsample(N,7);
			postSynNeuronIDs=randsample(N,divNum);%WITHOUT replacement (no double innervations)
			%postSynNeuronIDs=randsample(N,ceil(N/2));
			%postSynNeuronIDs=randsample(N,ceil(N));

			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%functionally, dendritic contacts can be categorized into 7 groups of input for a cell
			%assuming coactivity corresponds to similar dendritic output position targeting (perhaps developmentally)
			%assume a lot of such groups with various interconnections corresponding to sequence permutations
			%thus each neuron/subpopulation represents one random permutation of the active inputs (1-7)
			%***functional connectivity assumption (inactive connections inconsequential)*** 
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%presynNeuronIDs=allSequences(randi(5040),:);
			%presynNeuronIDs=allSequences(randi(numDendrites),:);
			%presynNeuronIDs=allSequences(randi(N),:);
			%seqConnMatrix(postSynID,:)=presynNeuronIDs;
			%seqConnMatrix(postSynID,presynNeuronIDs)=activitySeq;
			%randSeq=allSequences(randi(5040),:);
			%seqConnMatrix(postSynID,:)=randSeq; %assign a random sequence to the dendrite of this neuron
			
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%each neuron sends input to divNum neurons in population with random dendritic contact distance 
			%resulting in neurons with responses to many unique sequences in the other neurons
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			if(usePointConnectivity==0)
				seqConnMatrix(postSynNeuronIDs,preSynID)=randi(7,1,length(postSynNeuronIDs)); 
				
			else
				seqConnMatrix(postSynNeuronIDs,preSynID)=randi(1,1,length(postSynNeuronIDs)); 
			end
			%rMatrix=corr([activitySeq(:) presynNeuronIDs(:)], 'Type','Spearman');
			%rMatrix=corrcoef(activitySeq,presynNeuronIDs);
			%rMatrix=corrcoef(activitySeq,randSeq);
			%postSynTemplateR(postSynID)=rMatrix(2,1);
		end

		%use correlated connectivity assumption (e.g. Science dendritic sequence ripple paper, 2020)
		%{
		for i=2:(N-1)
		    for j=2:(N-1)
			p=rand(1);
			if(p<0.3 && seqConnMatrix(i,j)>1 && seqConnMatrix(i,j)<7 )
				seqConnMatrix(i,j-1)=seqConnMatrix(i,j)-1;
				seqConnMatrix(i,j+1)=seqConnMatrix(i,j)+1;
			end
		    end
		end
		%}

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%add somewhat local symmetric connectivity constraint
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		if(useLocalConnectivity)
			for i=1:N
			    for j=1:N
				if(i>=localConnectivityRadius && i<=N-localConnectivityRadius && abs(i-j)>localConnectivityRadius)
				%if(abs(i-j)>localConnectivityRadius && abs(i-(N-j+1))>localConnectivityRadius)
					seqConnMatrix(i,j)=0;
				end

				if(i<localConnectivityRadius && (j-i) < (N-localConnectivityRadius) && abs(i-j)>localConnectivityRadius)
					seqConnMatrix(i,j)=0;
				end

				if(i>N-localConnectivityRadius && (i-j) < (N-localConnectivityRadius) && abs(i-j)>localConnectivityRadius)
					seqConnMatrix(i,j)=0;
				end
			    end
			end
		end
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%remove self connections
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		if(noSelfConnections)
			for i=1:N
				seqConnMatrix(i,i)=0;
			end
		end
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%analyze expected dynamic range of sequence response by coverage of permutations
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



		%activeNeuronSeq=1:7;
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%sort connectivity matrix rows (post-synaptic neuron input vectors) 
		%by how much they are activated by initial input sequence
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%activeNeuronSeq=(N-6):N;
		%{
		activeNeuronSeq=1:7;
		for postSynID=1:N
			%rMatrix=corrcoef(activitySeq,randSeq);
			%rMatrix=corrcoef(activeNeuronSeq,seqConnMatrix(postSynID,activeNeuronSeq));
			
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%the timing of a subpopulation is considered its mode, representing gamma synchronization by I neurons 
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			
			substr1=sprintf('%d',1:7);
			substr2=sprintf('%d',seqConnMatrix(postSynID,activeNeuronSeq));
			%postSynTemplateR(postSynID)=rMatrix(2,1);
			[D dist str] = LCS(substr1,substr2);
			lcsLength=dist;	
			postSynTemplateR(postSynID)=lcsLength;
		end
		postSynTemplateR

		[~,sortedRidxes]=sort(postSynTemplateR(:));
		%seqConnMatrix=seqConnMatrix(sortedRidxes,:);
		%seqConnMatrix=seqConnMatrix(sortedRidxes,:);
		%}

		%figure
		imagesc(seqConnMatrix)
		%omarPcolor(1:N,1:numCompartmentsPerNeuron,(seqConnMatrix/max(seqConnMatrix(:)))')
		%omarPcolor(1:N,1:N,(seqConnMatrix/max(seqConnMatrix(:)))')
		%omarPcolor(1:N,1:N,(seqConnMatrix/max(seqConnMatrix(:))))
		ylabel('Post-synaptic Neuron No.')
		xlabel('Pre-synaptic Neuron No.')
		daspect([1 1 1])
		%ylabel('Dendritic Compartment No.')
		cb=colorbar
		ylabel(cb,'Synapse Dendritic Distance')
		%caxis([1 7])
		%daspect([1 1 1])
		%title('Random neural sequence templates (sorted by order)')
		title(sprintf('Population connectivity, %s',removeUnderscores(connectivityStr)))
		setFigFontTo(16)
		saveas(gcf,sprintf('../figures/connectivityMatrix_%s.tif',connectivityStr))

		%%xlim([1470 1550])
		%%ylim([1470 1550]+500)
		%ylim([2000 2200])
		%%title(sprintf('Subpopulation connectivity, %s (zoom)',removeUnderscores(connectivityStr)))
		%%saveas(gcf,sprintf('connectivityMatrixExZoom_%s.tif',connectivityStr))
		%fds
		%monitorSynapsesFig=figure;
		%monitorNetworkFig=figure;

				longSpikingNeuronIDs=[];

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%update activity state to next thetacycle 
		%based on 1) phase precession and 2) dendritic template matches
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%for cycleID=1:(numCycles-1)
		cycleID=0;
		hasMaxedOutBefore=false(N,1);
		visitCount=zeros(N,1);
		allCellsReached=false;
		fracs=[];

		isCurrentlyRefractory=false(N,numCycles);
		%refracPeriodNumCycles=20;
		refracPeriodNumCycles=0;
		%refracPeriodNumCycles=10;%1 sec refracrion after maxing out
		%refracPeriodNumCycles=1000;
		%refracPeriodNumCycles=100;
		
		%refracPeriodNumCycles=150;
		%refracPeriodNumCycles=200;
		%refracPeriodNumCycles=300;
		%refracPeriodNumCycles=500;
		%refracPeriodNumCycles=400;
		%refracPeriodNumCycles=5000;

		while(~allCellsReached)
			cycleID=cycleID+1;

			currSynapticActivity=zeros(N,N);

			%a(:,cycleID)=a(:,cycleID)/norm(a(:,cycleID));
			%postSynA(:,cycleID)=postSynA(:,cycleID)/norm(postSynA(:,cycleID));
				%%%%%%%%%%%%%%%%%%%%%%%%%
				%adaptation representation
				%%%%%%%%%%%%%%%%%%%%%%%%%
				if(cycleID>=5)
					%longSpikingNeuronIDs=[longSpikingNeuronIDs ;find(a(:,cycleID)>spikeThreshold & a(:,cycleID-1)>spikeThreshold & a(:,cycleID-2)>spikeThreshold & a(:,cycleID-3)>spikeThreshold)];
					longSpikingNeuronIDs=[longSpikingNeuronIDs ;find(a(:,cycleID)>spikeThreshold &...
							 a(:,cycleID-1)>spikeThreshold & a(:,cycleID-2)>spikeThreshold &...
							 a(:,cycleID-3)>spikeThreshold & a(:,cycleID-4)>spikeThreshold)];
					
					%a(longSpikingNeuronIDs,cycleID)=0;
					%a(longSpikingNeuronIDs,cycleID+1)=0;
				end
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%phase precession update; based on spiking threshold, whereas recruitment is based on top n winners
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%precessingCellIDs=(a(:,cycleID)>1);
			%precessingCellIDs=(a(:,cycleID)>minActivity); 
			%precessingCellIDs=(a(:,cycleID)>spikeThreshold); 
			precessingCellIDs=(a(:,cycleID)>spikeThreshold); 
			%a(:,cycleID+1)=a(:,cycleID)+activityPrecessionPerCycle;
			a(precessingCellIDs,cycleID+1)=a(precessingCellIDs,cycleID)+activityPrecessionPerCycle;
			%a(precessingCellIDs,cycleID+1)=a(precessingCellIDs,cycleID)+activityPrecessionPerCycle/2;

			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%how to permit currently non-active cell to join bump? 
			%not through phase precession (only for active cells including turn off)
			%but through recurrent network excitation selecting cells of some correlation
			%to current sequence; more cells on edge of bump than towards inside
			%only active cells precessing but all cells receiving functional network input 
			%which is WEAK relative to phase precession related input! 
			%(show this falls apart if network based movement is too strong?) 
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%saturation reset and normalize vector (simulating timing code feature)
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%a(:,cycleID+1)=a(:,cycleID+1)/norm(a(:,cycleID+1));
			%a(a>maxActivity)=0;
			maxedOutIDs=a(:,cycleID)>maxActivity;
			a(maxedOutIDs,cycleID)=0;
			a(maxedOutIDs,cycleID+1)=0;

			isCurrentlyRefractory(maxedOutIDs,cycleID:(cycleID+refracPeriodNumCycles))=true;
			hasMaxedOutBefore(maxedOutIDs)=true;
			visitCount(maxedOutIDs)=visitCount(maxedOutIDs)+1;
			fracs=[fracs length(find(hasMaxedOutBefore))/length(hasMaxedOutBefore)];
		
			%{	
			figure(cellCoverageFig)
			subplot(1,2,1)
			plot(fracs)
			xlabel('Cycle no.')
			title('fraction cells reached so far')
			subplot(1,2,2)
			bar(visitCount)
			xlabel('Cell no.')
			title('visit count per cell')
			%}
					

			%bar(hasMaxedOutBefore)
			if(isempty(find(~hasMaxedOutBefore)))
				allCellsReached=true;
			end
			%a(:,cycleID)=mod(a(:,cycleID),maxActivity+activityPrecessionPerCycle);
			%a=mod(a,maxActivity+2*activityPrecessionPerCycle);

			%biophysical threshold mechanism
			spikingNeuronIDs=find(a(:,cycleID)>spikeThreshold);

			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%winner takes all circuit spiking mechanism; most activated population sends dendritic potentials (least active silenced)
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			[sortedA,sortedIDs]=sort(a(:,cycleID));
			%spikingNeuronIDs=sortedIDs((end-assemblyN):end);	
			a(sortedIDs(1:(end-assemblyN-1)),:)=0;
			%a(sortedIDs(1:(end-assemblyN-1)),:)=spikeThreshold;

			

			disp(sprintf('%d neurons spiking...', length(spikingNeuronIDs)))
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%discretize activity levels for comparison to dendritic distances (~gamma synchrnoization)
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			discreteA=NaN(size(a(:,cycleID)));

			for d=1:8
				discreteA(a(:,cycleID)>=(d-1)*(maxActivity/8) & a(:,cycleID)<(d)*(maxActivity/8))=d-1;
			end	
			%figure; plot(a(:,cycleID)*7/max(a(:,cycleID))); hold on; plot(discreteA(:))

			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%sequence match update amongst neurons downstream of winner take all subpopulation
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%for dendriteNum=1:numDendrites
			for neuronID=1:N
				neuralActivityVectorMovie(neuronID,cycleID)=a(neuronID,cycleID);
				%if(~isempty(longSpikingNeuronIDs) && min(abs(longSpikingNeuronIDs-neuronID))==0)
				%	continue; % don't update adapting neurons (e.g.K channel blocks synaptic influence)
				%end
			%for neuronID=1:postSynN
				%for each neuron population unit, add in the similarity of template to actual  inputs for next cycle
				%templateSequence=seqConnMatrix(dendriteNum,:);
				%templateSequence=seqConnMatrix(neuronID,:);
				%observedSequence=a(:,cycleID);
				%observed sequence is actual activity sequence in average of each 7th of total pop 
				%{
				for i=1:7
					%startNeuron=(i-1)*floor(N/7)+1;
					%endNeuron=min(N,startNeuron+floor(N/7));
					startNeuron=N-(7-i+1)*assemblyN;
					endNeuron=min(N,startNeuron+assemblyN);
					observedSequence(i)=mean(a(startNeuron:endNeuron,cycleID));
				end
				%}
				%templateMatch=dot(templateSequence/norm(templateSequence),observedSequence/norm(observedSequence));

				presynNeuronTimings=discreteA(spikingNeuronIDs);
				presynNeuronTimings(isnan(presynNeuronTimings))=0;

				%fds
				%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
				%compute similarity between current activity vector and this neuron's template as count of timing/distance matches
				%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				templateMatch=0;	
				if(length(spikingNeuronIDs)>0)	
					%substr1=sprintf('%d',1:7);
					substr1=sprintf('%d',presynNeuronTimings);
					substr2=sprintf('%d',seqConnMatrix(neuronID,spikingNeuronIDs));
					%postSynTemplateR(postSynID)=rMatrix(2,1);
					%[D dist str] = LCS(substr1,substr2);
					%lcsLength=dist;
					%templateMatch=lcsLength
					%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
					%binary/locally-saturating synapses - ea timing-dendritic distance match only contributes 1 psp per dist
					%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
					%templateMatch=sum(substr1==substr2 & substr1~='0' & substr2~='0');
					templateMatch=0;
					for d=1:7
						currDistChar=num2str(d);
						numMatchesCurrDist=sum(substr1==substr2 & substr1==currDistChar & substr2==currDistChar);
						if(numMatchesCurrDist>0)
							hasMatchForCurrDist=1;
						else
							hasMatchForCurrDist=0;
						end

						templateMatch=templateMatch+hasMatchForCurrDist;
					end
				end
				%for s=1:length(spikingNeuronIDs)
				%	currSynapticActivity(neuronID,spikingNeuronIDs(s))=seqConnMatrix(neuronID,spikingNeuronIDs(s));
				%end
			
				%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				%record active synapses in this cycle and their weights
				%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				currSynapticActivity(neuronID,spikingNeuronIDs(substr1==substr2))=seqConnMatrix(neuronID,spikingNeuronIDs(substr1==substr2));
				%neuronID=ceil((dendriteNum/numDendrites) *N);
				%a(neuronID,cycleID+1)=a(neuronID,cycleID)+templateMatch*networkExcStrength;
				%a(neuronID,cycleID+1)=a(neuronID,cycleID)+activityPrecessionPerCycle*(templateMatch/7);
				%a(neuronID,cycleID+1)=a(neuronID,cycleID)+networkExcStrength*(templateMatch/7);
				%if(templateMatch>4)
		
				%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				%update activity in next cycle based on template match based dendritic drive in current cycle 
				%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				currDendriticDrive=0;
				%refraction prevents further input via inhibition
				refractFactor=1;
				if(isCurrentlyRefractory(neuronID,cycleID))
					%templateMatch=templateMatch/2;			
					%templateMatch=0;
					refractFactor=0.2;
				end

				if(sequenceMatchDriven)
					%sequence match-based dendritic integration	
					if(templateMatch>=2)
					%if(templateMatch>=4)
						%in order to re-create sequence/gradient in next cycle,need large dynamic range of response
						currDendriticDrive=refractFactor*networkExcStrength*(templateMatch/7);
						
						%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
						%how many neurons are expected to be recruited per cycle? Sustained?
						%relative importance of initial recruitment vs cross-cycle sustainment? (later seems more)
						%rate of phase precession per cycle vs magnitude of cell recruitment per cycle?
						%how are these coupled?
						%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
						a(neuronID,cycleID+1)=a(neuronID,cycleID)+currDendriticDrive;
						%a(neuronID,cycleID+1)=a(neuronID,cycleID)+networkExcStrength*(templateMatch/2);
						%a(neuronID,cycleID+1)=a(neuronID,cycleID)+networkExcStrength*(templateMatch);
					end
				else		
					%linear dendritic integration control
					a(neuronID,cycleID+1)=a(neuronID,cycleID)+networkExcStrength/100*sum(seqConnMatrix(neuronID,spikingNeuronIDs));		
				end

				
				%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				%record total synaptic drive to this neuron in current cycle
				%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				dendriticActivationOverTime(neuronID,cycleID)=currDendriticDrive;

				%a(neuronID,cycleID+1)=a(neuronID,cycleID)+seqConnMatrix(neuronID,:) *networkExcStrength;
				%a(:,cycleID+1)=a(:,cycleID)+seqConnMatrix*a(:,cycleID)*networkExcStrength;
				%postSynA(:,cycleID+1)=postSynA(:,cycleID)+seqConnMatrix*a(:,cycleID);
				%postSynA(:,cycleID+1)=postSynA(:,cycleID)+seqConnMatrix*a(:,cycleID);
			end %neuron id loop
			%a(:,cycleID+1)=a(:,cycleID+1)/norm(a(:,cycleID+1));
			%postSynA(:,cycleID+1)=postSynA(:,cycleID+1)/norm(postSynA(:,cycleID+1));

			%if(cycleID==maxCycleDisp)
			if(mod(cycleID,50)==0)
				figMov=figure;
				%omarPcolor(1:numCycles,1:N,neuralActivityVectorMovie,figMov)
				omarPcolor(1:cycleID,1:N,neuralActivityVectorMovie(:,1:cycleID),figMov)
				colormap(parula)
				cb=colorbar
				drawnow
				xlabel('Cycle num')
				ylabel('Neuron no.')
				ylabel(cb,'Activity level')
				%saveas(gcf,sprintf('circuitActivityOverTime_%s.tif',connectivityStr))	
				%fds
			end
			
			%figure(monitorSynapsesFig)
			%figure(monitorNetworkFig)
			%subplot(2,1,1)

			%for n=1:N
			%	currSynapticActivity(n,n)=NaN;
			%end

			%functionalConnectivityMovie(:,:,cycleID)=currSynapticActivity;		

			%{
			imagesc2(currSynapticActivity)
			ylabel('Post-synaptic Neuron No.')
			xlabel('Pre-synaptic Neuron No.')
			%daspect([1 1 1])
			%ylabel('Dendritic Compartment No.')
			cbA=colorbar
			ylabel(cbA,'Synaptic activity')
			%caxis([1 7])
			%daspect([1 1 1])
			%title('Random neural sequence templates (sorted by order)')
			title(sprintf('Functional synapses, %s',removeUnderscores(connectivityStr)))

			caxis([0 7])
			subplot(2,1,2)
			
			%subplot(1,2,1)
			%plot(a(1:N,cycleID),'o','Color',cycleCMap(cycleID,:))
			%hold on
			plot(a(1:N,cycleID),'-','Color',cycleCMap(cycleID,:))
			%hold on

			xlabel('Neuron No.')
			ylabel('Activity level')
			title({sprintf('Sequence-length driven network across theta cycles, %s',removeUnderscores(connectivityStr)),sprintf('Cycle %d',cycleID)})
			maxFigManual2d(1,1,16)
			%setFigFontTo(16)
			cb2=colorbar
			ylabel(cb2,'Time (a.u.)')

			ylim([0 maxActivity*1.1])	
			%subplot(1,2,2)
			%plot(postSynA(:,cycleID),'Color',cycleCMap(cycleID,:))
			%hold on
			%maxFig
			%drawnow 
			%pause(0.15)
			%}
			
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%break conditions
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

			%hit enough cells to mostly determine ordering
			if(length(find(hasMaxedOutBefore))>= minNumVisitedCells)
				numCycles=cycleID;
				break
			end
			%no more potential for activity revival
			if(max(a(:,cycleID))==0)
				numCycles=cycleID;
				break
			end
		end %theta cycle loop

		xValues=1:N;
		A=neuralActivityVectorMovie(:,1:cycleID);
		%A=A.*(A>0.8);
		columnDescrip='placeCellActivityVector'

		xUnitStr='Cell no.';
		yUnitStr='Activity';
		iterationUnitStr='Cycle num';

		%svdAccessRoutine(xValues,A,columnDescrip,xUnitStr,iterationUnitStr,yUnitStr)

		%nnmfAccessRoutine(xValues,A,columnDescrip,xUnitStr,iterationUnitStr,yUnitStr)
		
		if(detectAssemblies)
			figComps=figure(111)
			figTimes=figure(112)
			[W,H]=seqNMFaccessRoutine(xValues,A,columnDescrip,xUnitStr,iterationUnitStr,yUnitStr,figComps,figTimes)

			figure(figComps)
			maxFigManual2d(1/3,2.8/3,4)
			
			figure(figTimes)
			maxFigManual2d(1/3,2.8/3,4)

			saveas(figComps,sprintf('../figures/detectedPlaceAssemblies_%s_PhasePrecessionSpeed_%.2f.tif',connectivityStr,speedFact))	
			saveas(figTimes,sprintf('../figures/detectedPlaceAssemblyActivationCycles_%s_PhasePrecessionSpeed_%.2f.tif',connectivityStr,speedFact))	
			save(filePath,'rngVal','frozenUniformNoiseConditions','W','H','A','speedFact','connectivityStr','functionalConnectivityMovie','dendriticActivationOverTime','numCycles','-v7.3')
		else
			if(plotDendriticDrive)
				save(filePath,'rngVal','frozenUniformNoiseConditions','A','speedFact','connectivityStr','functionalConnectivityMovie','dendriticActivationOverTime','numCycles','noiseConditionID','-v7.3')
			else
				save(filePath,'rngVal','frozenUniformNoiseConditions','A','speedFact','connectivityStr','numCycles')
			end
		end
		%trajectory is distance from each component at each time
		%distance representation has to do with mutual similarity between columns of W
		%numFactorsUsed=17;
		
	else
		load(filePath)
	end
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%assign each neuron to place assembly chosen by activity flow
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	figure(cellCoverageFig)
	subplot(1,2,1)
	plot(fracs)
	xlabel('Cycle no.')
	title('fraction cells reached so far')
	subplot(1,2,2)
	bar(visitCount)
	xlabel('Cell no.')
	title('visit count per cell')
	saveas(gcf,sprintf('../figures/cellCoverageByActivityStream_%s_PhasePrecessionSpeed_%.2f.tif',connectivityStr,speedFact))

	activityThresh=0.1;

	numNeurons=size(A,1);
	neuronCompMap=zeros(numNeurons,1);
	neuronCOMmap=zeros(numNeurons,1);

	for neuronID=1:numNeurons
		%neuronCOMmap(neuronID)=getCenterOfMass(A(neuronID,:));
		%neuronCOMmap(neuronID)=neuronID;
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%only consider times of reaching top of field
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		highVals=A(neuronID,:);
		highVals(highVals<(maxActivity)*0.9)==0;
		[~,maxID]=max(highVals);
		%[~,maxID]=max(A(neuronID,:));
		neuronCOMmap(neuronID)=maxID;
	end

	for neuronID=1:numNeurons
		%[~,maxID]=max(A(neuronID,:));
		%neuronCOMmap(neuronID)=maxID;
	end

	if(detectAssemblies)
		Wnorm=normalizeCols(W);
		for i=1:size(Wnorm,2)
			currCompNeuronIDs=find(Wnorm(:,i)>activityThresh);
			neuronCompMap(currCompNeuronIDs)=i;
		end
		[~,sortByAllegianceThenCOM]=sortrows([neuronCompMap(:) neuronCOMmap(:)]);
	else
		[~,sortByAllegianceThenCOM]=sortrows( neuronCOMmap(:) );
	end

	%[~,sortByAllegiance]=sort(neuronCompMap,'descend');
	%[~,sortByAllegianceThenCOM]=sortrows(-[neuronCompMap(:) neuronCOMmap(:)]);

	%figure; imagesc(fliplr(A(sortByAllegiance,:)))
	%[~,sortByAllegiance]=sort(neuronCompMap);
	%figure; imagesc(A(sortByAllegiance,:))
	figH=figure;

	functionalDI=getDI(sortByAllegianceThenCOM);
	maxDI=getDI(minNumVisitedCells:-1:1);
	normDI=functionalDI/maxDI;

	
	normDIall(noiseConditionID)=normDI;
	rawDIall(noiseConditionID)=functionalDI;

	%save(filePath,'functionalDI','maxDI','normDI','-append')
	
	save(filePath,'rawDIall','maxDI','normDIall','-append')

	%omarPcolor(1:size(A,2),1:size(A,1),A(sortByAllegiance,:),figH)
	if(noiseConditionID==1)
		omarPcolor(1:size(A,2),1:size(A,1),flipud(A(sortByAllegianceThenCOM,:)),figH)
		ylabel('Neuron no. (sorted by component activation time)')
		xlabel('Cycle number')
		%ylim([1 250])
		%ylim([1 275])
		ylim([1 N])
		title(sprintf('Temporal rate of phase precession: %.2f',speedFact))
		saveas(gcf,sprintf('../figures/circuitActivityOverTime_%s_PhasePrecessionSpeed_%.2f.tif',connectivityStr,speedFact))	
	end

	%resortedConnMov=functionalConnectivityMovie(sortByAllegianceThenCOM,sortByAllegianceThenCOM,:);
	%resortedConnMov=flipud(functionalConnectivityMovie(sortByAllegianceThenCOM,sortByAllegianceThenCOM,:));

	if(saveActiveSynapseMovie)
		resortedConnMov=functionalConnectivityMovie(sortByAllegianceThenCOM,sortByAllegianceThenCOM,end:-1:1);
		movFilePath=sprintf('../figures/synapticActivityOverTime_%s_PhasePrecessionSpeed_%.2f.avi',connectivityStr,speedFact);
		frameRate=10;
		colormapName=parula
		matrix2MovieTibin(movFilePath,resortedConnMov,frameRate,colormapName)
	end


	if(plotDendriticDrive)
		startCellID=1;
		endCellID=N;
		resortedDendriticDrive=dendriticActivationOverTime(sortByAllegianceThenCOM,:);
		figure
		smoothSizeCycleIDs=5;
		%smoothSizeCellIDs=5;
		%smoothSizeCellIDs=13;
		smoothSizeCellIDs=13;
		%smoothSizeCycleIDs=4;
		%smoothSizeCycleIDs=6;
		%smoothSizeCycleIDs=2;
		frameCount=1;
		for cycleID=(smoothSizeCycleIDs+1):(numCycles-smoothSizeCycleIDs)
			%activityDistOverTime(:,cycleID)=sum(resortedConnMov(:,:,cycleID),1);
			%dendriticActivationOverTime(:,cycleID)=sum(resortedConnMov(:,:,cycleID),2);
			hold on
			%if(cycleID>1)
			subplot(2,1,1)
				%plot(A(sortByAllegianceThenCOM,(cycleID-smoothSizeCycleIDs):(cycleID+smoothSizeCycleIDs)),'k')
			%end

			currA=A(sortByAllegianceThenCOM,cycleID-1);
			currA=movmean(currA,smoothSizeCellIDs);	
			%[currCOM,~]=getCenterOfMass(currA(startCellID:endCellID));
			[~,currCOM]=max(currA(startCellID:endCellID));
			currCOM=currCOM+startCellID-1;
			
			plot(currA,'k')
			
			currDist=mean(resortedDendriticDrive(:,(cycleID-smoothSizeCycleIDs):(cycleID+smoothSizeCycleIDs)),2);

			currDist=movmean(currDist,smoothSizeCellIDs);	
			%currDist=movmedian(currDist,smoothSizeCellIDs);	

			%[currCOMdend,~]=getCenterOfMass(currDist(startCellID:endCellID));
			[~,currCOMdend]=max(currDist(startCellID:endCellID));
			currCOMdend=currCOMdend+startCellID-1;

			%ylim([0 1.7])
			ylim([0 1])
			xlim([startCellID endCellID])	
			hold on

			plot([currCOM currCOM],[0 1.7],'k--','LineWidth',2)
			plot([currCOMdend currCOMdend],[0 1.7],'r--','LineWidth',2)
			xlabel('Cell no.')
			ylabel('Activity level (~theta phase)')
			legend('activity')
			subplot(2,1,2)

			plot(currDist,'r')
			hold on
			%yMax=0.075;
			%yMax=0.3;
			%yMax=0.15;
			%yMax=0.2;
			try
				yMax=networkExcStrength/5;
			catch
				yMax=0.1;
			end
			plot([currCOM currCOM],[0 yMax],'k--','LineWidth',2)
			plot([currCOMdend currCOMdend],[0 yMax],'r--','LineWidth',2)
			
			xlabel('Cell no.')
			ylabel('dendritic drive')
			%ylim([0 0.15])
			ylim([0 yMax])
			xlim([startCellID endCellID])	
			title(sprintf('cycle num %d',cycleID))
			legend('dendritic activation')
			%pause(0.1)
			hold off
			maxFig
			setFigFontTo(16)	
			drawnow
			F(frameCount)=getframe(gcf)
			frameCount=frameCount+1;
			clf;
		end
		movFilePath=sprintf('../figures/synapticDriveVsSpikingActivity_%s_PhasePrecessionSpeed_%.2f.avi',connectivityStr,speedFact);
		frameArray2Movie(movFilePath,F)
	end
	
	disp(sprintf('finished computing for noise condition id %d.',noiseConditionID))
	toc


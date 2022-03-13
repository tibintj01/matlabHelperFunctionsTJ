classdef Cells < handle & matlab.mixin.Copyable %create object by reference
	%encapsulate data and actions of cells to keep my interface and implementation details separate
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%public cells properties
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	properties(Constant)
		justSpikingConductances=1
		%BLOCK_OUTPUT_SPIKING=1
		BLOCK_OUTPUT_SPIKING=0
		
		%ADJ_LIN_DELAY_BASELINE=40+190-285.714
		%ADJ_LIN_DELAY_BASELINE=-56+20
		ADJ_LIN_DELAY_BASELINE=-56+20+63
		ADJ_LIN_DELAY_SLOPE=5;
		ADJ_LOG_DELAY_BASELINE=20 %with ADJ_LOG_DELAY_SLOPE=1
		%ADJ_LOG_DELAY_BASELINE=20+150-72 %with ADJ_LOG_DELAY_SLOPE=2.5
		%ADJ_LOG_DELAY_SLOPE=2.5;
		ADJ_LOG_DELAY_SLOPE=1;
		includeKS=1;
		NUM_CELLS_L2=1;
		%INTEGRATOR_gL=0.005;
		%INTEGRATOR_gL=0.1;
		INTEGRATOR_gL=0.01;
		SEQ_DET_gL=0.1; %CA1 10msec time constant
		%SEQ_DET_gL=0.2; %CA1 10msec time constant
	end
	properties
		%voltage and time dependent gating variable matrices
		%network of phase precessing cells is both sensor and effector of internally-generated neural movement
		%synapses give salience but spatial information is encoded by theta phase (perhaps through inhibitory control)
		v

		injCurrMatrix
		inhThetaInputArray

		numCellsPerPlace
		numPlaces
		numSteps

		internalConnObj
		externalInputObj
		feedforwardConnObj
	
		dt

		spikeTimes
		delayedSpikeTimes
		doubleDelayedSpikeTimes
		spikeCellCoords
		gsyn
		
		IsynRecord
		EsynRecord
		
		esyn_I
		esyn_E

		%circuitRawOutput
		%sensoryInput
		%speedInput
		%placeInput
		%lfpInput	
	end
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%cells conductance parameters
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	properties(Access=public)
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%passive (linear) conductance parameters
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%specific membrane resistance 30,000 ohm-cm^2 (Gm=0.033mS/cm^2)
		% i.e. leak current
		%gl=0.033;        %mS/cm^2
		%gl=0.033*8;        %mS/cm^2
		%gl=0.033*5;        %mS/cm^2
		%gl=0.033*6;        %mS/cm^2
		%gl=SCAN_PARAM2;        %mS/cm^2
		gl=0.0333333333333;        %mS/cm^2

		el=-60;         %mV "normally -70 mV but adjusted to keep resting membrane potential near -66 mV." (Leung, 2011)
		cm=1.0; %uF/cm^2		

	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%public methods
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods
		%constructor
		function thisObj=Cells(simProps)
			%set properties of this cells object
			if (nargin == 1)
				thisObj.setSimSpecificProperties(simProps);

				%if(isfield(simProps,'bias'))
				%	thisObj.setCellProps(simProps.extEnvObj,simProps.bias);
				%else
					thisObj.setCellProps(simProps.extEnvObj);
				%end

				thisObj.setIntrinsicsMatrix();
				
				%thisObj.injCurrMatrix=thisObj.currInjectorArray.getFloatMatrix();
				thisObj.injCurrMatrix=thisObj.externalInputObj.getFloatMatrix();
				
				thisObj.spikeTimes=[];
				thisObj.delayedSpikeTimes=[];
				thisObj.doubleDelayedSpikeTimes=[];
				thisObj.spikeCellCoords=[];
				
				thisObj.gsyn=zeros(thisObj.numCellsPerPlace,thisObj.numPlaces,thisObj.numSteps);
				
				thisObj.esyn_E=thisObj.internalConnObj.esyn_E;
				thisObj.esyn_I=thisObj.internalConnObj.esyn_I;

			%elseif(nargin==2 && strcmp(initStr,'backbone'))
			%	thisObj.setSimSpecificProperties(simConfigObj);
                        %        thisObj.setCellProps();
			end
		end


		%stepping through time
		function go(thisObj)
			disp('solving cell diff equations....')
				thisObj.injCurrMatrix=thisObj.externalInputObj.getFloatMatrix();
			tic
			thisObj.letItRip()
			toc
		end
		
		
		function idStr=getCellIDstr(thisObj,r,placeIdx)
			itonic=thisObj.injCurrMatrix(r,placeIdx);
			idStr=sprintf('i_{tonic}=%.2f',gnap,gks,itonic);
		end
	end
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%private methods
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods(Access= protected)
		function letItRip(thisObj)
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%load variables		
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			numSteps=thisObj.numSteps;
			numPlaces=thisObj.numPlaces;
			numCellsPerPlace=thisObj.numCellsPerPlace;
			v=thisObj.v;
			
			nks=thisObj.nks;
			gl=thisObj.gl;
			el=thisObj.el;
			
			ek=thisObj.ek;
			
			includeKS=Cells.includeKS;	
			%make sure these are updated based on gbar and gsigma
			thisObj.setIntrinsicsMatrix();
			gksMatrix=thisObj.gksMatrix;		
			
			injCurrMatrix=thisObj.injCurrMatrix;
			dt=thisObj.dt;
			cm=thisObj.cm;

			gInhThetaMatrix=thisObj.inhThetaInputArray.conductanceTimeSeries;
			
			thetaPhaseSeries=thisObj.inhThetaInputArray.getPhaseOverTime(1,1);

			distToPhaseFact=360/(ExternalEnvironment.PLACE_INPUT_WIDTH/ExternalEnvironment.NUM_CYCLES_PER_FIELD);
			
			spikeTimes=thisObj.spikeTimes;
			delayedSpikeTimes=thisObj.delayedSpikeTimes;
			doubleDelayedSpikeTimes=thisObj.doubleDelayedSpikeTimes;
			spikeCellCoords=thisObj.spikeCellCoords;

			gsyn=thisObj.gsyn;
			
			esyn_I=thisObj.esyn_I;
			esyn_E=thisObj.esyn_E;

			tausyn=thisObj.internalConnObj.tausyn;
			connectivityMatrix=thisObj.internalConnObj.connectivityMatrix;

			feedfwdGmatrix=thisObj.feedforwardConnObj.connectivityMatrix;
                        dendriticDelayTemplateMatrix=thisObj.feedforwardConnObj.dendriticDelayTemplateMatrix;

			startCouplingTime=thisObj.internalConnObj.startCouplingTime;

			iksRecord=NaN(size(v));

			IsynRecord=NaN(size(vL2));
			EsynRecord=NaN(size(vL2));
			
			normFactor=DelayObject.NORM_FACTOR;
			convFactor=DelayObject.CONV_FACTOR;
			baselineDelay=DelayObject.BASELINE_DELAY;	
			imax=DelayObject.IMAX;
			imin=DelayObject.IMIN;


			minPrecession=thisObj.feedforwardConnObj.minDelay;
			if(FeedForwardConnectivity.USE_LINEAR_DELAYS)
				maxPrecession=thisObj.feedforwardConnObj.maxDelay;	
				%linearPrecessionSlope=(maxPrecession-minPrecession)/(imax-imin);
				%linearPrecessionSlope=2*(maxPrecession-minPrecession)/(imax-imin); %real phase precession ends at 180 degrees-> double slope?
				%linearPrecessionSlope=3*(maxPrecession-minPrecession)/(imax-imin); %real phase precession ends at 180 degrees-> double slope?
				%linearPrecessionSlope=5*(maxPrecession-minPrecession)/(imax-imin); %real phase precession ends at 180 degrees-> double slope?
				linearPrecessionSlope=Cells.ADJ_LIN_DELAY_SLOPE*(maxPrecession-minPrecession)/(imax-imin); %real phase precession ends at 180 degrees-> double slope?
			end
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%step through time
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			for step=1:numSteps-1
				if(mod(step,100000)==1)
					disp(sprintf('INTEGRATING!!!!! t=%d msec', round(dt*step)))
					%disp(sprintf('INTEGRATING!!!!! v=%.5f',v(1,1,step)))
				end
				for placeIdx=1:numPlaces
					for cellRow=1:numCellsPerPlace
						vSpecific=v(cellRow,placeIdx,step);
						itonic=injCurrMatrix(cellRow,placeIdx,step);

						if(includeKS)
							if(step==1)
								nks(cellRow,placeIdx,step)=0.5; %%%%
							end
							nksSpecific=nks(cellRow,placeIdx,step);
						end
						%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
						%check for spike in  time step; resets V and sends a PSP into the next cycle at phase; spikes due to phase precession (timing in last cycle) 
						%vs spikes due to dendritic integration (intensity in this cycle) - couldn't these be in conflict? Inputs in this cycle vs delayed influence from last cycle?
						%Perhaps this implies conditions of self-consistency between actual timing of network input on next cycle and predicted timing by phase precession inhib network?
						%well, since predicted timing uses consistent speed info (via sequence response phase) and network timing is controlled by coordinated precession, these should match!
						%what is their relative contribution to the spike in the next cycle? Half and half? PSP (network) vs release of inhibition (phase precession)? 
						%let's start with half and half psp strengths but make the relative strengths a parameter to analyze its role later
						%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
						 
						%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
						%phase precession is for spatial updating (speed integration), dendritic integration is for sequence detection... spatial updating based on sequence match?
						%phase precession proportional to # inputs in order with memory across multiple theta cycles (to produce asymmetry; perhaps manifesting as faster oscillator)
						%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
						%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
						%a spike happens; 
						%-immediately send a psp to all post-synaptic partners with delays corresponding to contact distances; 
						%-also send a phase precession PSP to this cell in next cycle based on timing in this cycle;
						%only when internally predicted timing and network dendritic sequence timing match up will cell keep spiking;    
						%***set spike threshold at voltage corresponding to around "3-4 in order" response***
						%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
						if(step>1 && v(cellRow,placeIdx,step)>-30 && v(cellRow,placeIdx,step-1) <-30)
							spikeTimes=[spikeTimes; step*dt];
							if(numCellsPerPlace*numPlaces>=1)
								cellCoord=[cellRow placeIdx];
								spikeCellCoords=[spikeCellCoords; cellCoord];
								
								depConstant=1;
								numRecentSpikes=getNumSpikesInLastWind(spikeTimes,spikeCellCoords,FeedForwardConnectivity.SYN_DEP_WINDOW);
								 depConstant=(FeedForwardConnectivity.SYN_DEP_FACT)^(numRecentSpikes);

								%get synaptic conductance time course for all cells
								 %if this cell spikes, add a synaptic weight time course to all of
								 %its post-synaptic recipients' gsyn
								 %400 ms covers integrated timecourses without slowing down
								 %synEndStep=step+1+round(400/dt);
								 %4000 ms covers integrated timecourses with delay
								 synEndStep=step+1+round(4000/dt);
								 if(synEndStep>numSteps)
								     synEndStep=numSteps;
								 end
								 synCurrentIdxes=(step+1):synEndStep;

								 %count=0;
								%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                                %L1 post-synaptic conductance changes
                                                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                                for postSynCellIdx=1:numCellsPerPlace
                                                                         for postSynPlaceIdx=1:numPlaces
                                                                             weight=connectivityMatrix(cellRow,placeIdx,postSynCellIdx,postSynPlaceIdx);

                                                                             if(step*dt>=startCouplingTime && weight>0)
                                                                                gsyn(postSynCellIdx,postSynPlaceIdx,synCurrentIdxes)=squeeze(gsyn(postSynCellIdx,postSynPlaceIdx,synCurrentIdxes))...
                                                                                    +weight*exp(-dt*(synCurrentIdxes-(step+1))/tausyn).'; %instantaneous conductance jump with single additive exp decay..
                                                                                %count=count+1;
                                                                             end
                                                                        end
                                                                end
								     
                                                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                                %inhibitory delayed post-synaptic conductance changes
								%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
								     defaultPhaseSlope=baselineDelay/(imax-imin);
                                                                     
									%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
									%amount of delay due to phase precession is a function of current theta phase and current speed (through inhibition and excitation signals)
									%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
									currentThetaPhase=thetaPhaseSeries(step);
									currentSpeed=ExternalEnvironment.CONSTANT_RUN_SPEED;

									%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
									%how to ensure this ends up being a logarithmic or linear phase vs position graph? Have the correct phase derivatives wrt time and position
									%ie dPhase/dt=k*(currSpeed) --> temporal phase prec. slope m=k*(currSpeed) (k=distance to phase conversion)
									%cycle to cycle phase update: LOG: dPhase/dx=m/x--> phase2=phase1+m/phase1, LINEAR: dPhase/dx=m (on this time scale, dx=speed*dt) --> phase2=phase1+m
									%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
									%phasePrecessionDelay=Cells.ADJ_LOG_DELAY_SLOPE*(normFactor*log(convFactor*(imax-itonic))+minPrecession-Cells.ADJ_LOG_DELAY_BASELINE)+(defaultPhaseSlope*(itonic-imin)); %see DelayObject for values
									tempPhasePrecSlope=distToPhaseFact*(currentSpeed);
									
									phaseDeltaPerCycle=tempPhasePrecSlope*(frequencyDefault/1000);
									%log phase vs position relationship
									phasePrecessionDelay=phaseDeltaPerCycle/currentThetaPhase;

								     if(FeedForwardConnectivity.USE_LINEAR_DELAYS)
									%phasePrecessionDelay=-(linearPrecessionSlope*(itonic-imin)+minPrecession);	
									%phasePrecessionDelay=-(linearPrecessionSlope*(itonic-imin)+minPrecession)+(defaultPhaseSlope*(itonic-imin));
									%cancel out model precession and add in synthetic phenomenological precession	
									%phasePrecessionDelay=(linearPrecessionSlope*(imax-itonic)+minPrecession-Cells.ADJ_LIN_DELAY_BASELINE)+(defaultPhaseSlope*(itonic-imin));	
									%linear phase vs position relationship
								     	phasePrecessionDelay=phaseDeltaPerCycle;
								     end
									%phasePrecessionDelay=normFactor*log(convFactor*(imax-itonic))+baselineDelay; %see DelayObject for values
									%edge cases
                                                               		%if(~isreal(phasePrecessionDelay) || phasePrecessionDelay<0)
                                                               		if(~isreal(phasePrecessionDelay)) %often negative, as log argument becomes imaginary!!
										phasePrecessionDelay=0;
									end 
                                                        	    delayedSpikeTimes=[delayedSpikeTimes ;(step*dt + phasePrecessionDelay)];

							end%more than one cell condition
						end%spike condition

						%%%%%%%%%%%%%%%%%%
						%Il, leak current
						%%%%%%%%%%%%%%%%%%
						il=gl*(vSpecific-el);

						if(includeKS)
						      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
						      %I_KS gate voltage and time dependence
						      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
						      %gks=8;
							%{
							if(cell==1)
								gks=gks1;
							elseif(cell==2)
								gks=gks2;
							end
							%}
							gks=gksMatrix(cellRow,placeIdx);

							nks_Vt=-35;
							%nks_Vt=-45;
							%nks_Vt=-50;
							%nks_Vt=-55;
						      %nks_Vt=-10;
						      %nks_Gain=-10;
						      nks_Gain=-11;
						      nksInf=1/(1+exp((vSpecific-nks_Vt)/nks_Gain));


						       tau_nks=1/((exp((vSpecific-nks_Vt)/40)+exp(-(vSpecific-nks_Vt)/20) )/81);

							%tau_nks=tau_nks/2;

						      %slow, low threshold potassium current
						      iks=gks*nksSpecific*(vSpecific-ek);
						end


						%synaptic current
						%{
						if(numCells>1)
							isyn=gsyn(step,cell)*(vSpecific-esyn);
						else
							isyn=0;
						end
					
						isynExt=g_Inh(step)*(vSpecific-esynI) + g_Exc(step)*(vSpecific-esynE);
						%}

						isynIntE=gsyn(cellRow,placeIdx,step)*(vSpecific-esyn_E);
						isynExt=gInhThetaMatrix(cellRow,placeIdx,step)*(vSpecific-esyn_I);
					      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
					      %increment variables using euler's method of ODE integration
					      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
						iksRecord(cellRow,placeIdx,step)=iks;
						
						%vInc=double(dt*(-il-ina-ik-ika-ih-inap-iks-isynIntE-isynExt+itonic)/cm);
						vInc=double(dt*(-il-iks-isynIntE-isynExt+itonic)/cm);

					       %V=V+dV
					      v(cellRow,placeIdx,step+1)=vSpecific+vInc;

						if(includeKS)
							nksInc=double(dt*(nksInf-nksSpecific)/tau_nks);
						      nks(cellRow,placeIdx,step+1)=nksSpecific+nksInc;
						end
					end %loop over cells coding for current place
				end %loop over places
			end %loop over time steps
			
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%store raw output
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		 	thisObj.v=v;
		 	%thisObj.gsynL2=gsynL2;
		 	thisObj.gsyn=gsyn;
		 	
			%thisObj.l2EsynRecord=l2EsynRecord;
		 	%thisObj.l2IsynRecord=l2IsynRecord;
			thisObj.EsynRecord=EsynRecord;
		 	thisObj.IsynRecord=IsynRecord;
			

			thisObj.nks=nks;
			%thisObj.inaRecord=inaRecord;
			%thisObj.ikRecord=ikRecord;
			thisObj.iksRecord=iksRecord;
			thisObj.spikeTimes=spikeTimes;
			thisObj.delayedSpikeTimes=delayedSpikeTimes;
			thisObj.doubleDelayedSpikeTimes=doubleDelayedSpikeTimes;
			thisObj.spikeCellCoords=spikeCellCoords;
		end %letItRip function

		function setIntrinsicsMatrix(thisObj)
			thisObj.gnapMatrix=normrnd(thisObj.gnapBar,thisObj.gnapSigma,thisObj.numCellsPerPlace,thisObj.numPlaces);
		        thisObj.gksMatrix=normrnd(thisObj.gksBar,thisObj.gksSigma,thisObj.numCellsPerPlace,thisObj.numPlaces);
			thisObj.gnapMatrix(thisObj.gnapMatrix<0)=0;
			thisObj.gksMatrix(thisObj.gksMatrix<0)=0;
			
			%sort in order of increasing excitability
			thisObj.gnapMatrix=sort(thisObj.gnapMatrix,'descend');
			thisObj.gksMatrix=sort(thisObj.gksMatrix,'ascend');
		end	
		
		function setSimSpecificProperties(thisObj,simParams)
			thisObj.numCellsPerPlace=simParams.numCellsPerPlace;
			thisObj.numPlaces=simParams.numPlaces;
			thisObj.numSteps=simParams.numSteps;
			thisObj.dt=simParams.dt;
		end


		function setCellProps(thisObj,extEnvObj,bias)
			numCellsPerPlace=thisObj.numCellsPerPlace;
			numPlaces=thisObj.numPlaces;
			numSteps=thisObj.numSteps;
			dt=thisObj.dt;
			
			timeAxis=(dt:dt:(numSteps*dt)).'; %assumes simTime is divisible by dt
			
			simSpecificInfo.numCellsPerPlace=numCellsPerPlace;
			simSpecificInfo.numPlaces=numPlaces;
			simSpecificInfo.timeAxis=timeAxis;
			if(exist('bias','var'))
				simSpecificInfo.bias=bias;
			end

			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%cell settings
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%gksBar=0.2;
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%Oct 2 2019, TJ
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%gksBar=0.1;
			%gksBar=0.15;
			%gksBar=0.3;
			%gksBar=0.6;
			%gksBar=0.7;
			%gksBar=0.8;
			%gksBar=1.6;
			%gksBar=2.1;
			%gksSigma=0.05;
			
			%Nov 9 2019 - TJ, reduce variability
			gksBar=2;
			%gksSigma=0.01;
			%gksSigma=0.005;
			%gksSigma=0.001;
			gksSigma=0;
			
			%gksSigma=0.1;

			if(numCellsPerPlace*numPlaces==1)
				gksSigma=0;
			end

			%gnapMatrix=normrnd(gnapBar,gnapSigma,numCellsPerPlace,numPlaces);
			%gksMatrix=normrnd(gksBar,gksSigma,numCellsPerPlace,numPlaces);

			currInjectorArray=CurrentInjectors(simSpecificInfo,extEnvObj);
			%currInjectorArray.displayContent();
			%drawnow
			%injCurrMatrix=currInjectorArray.getFloatMatrix();


			inhThetaInputArray=ThetaPopInput(simSpecificInfo);
			%inhThetaInputArray.displayContent();


			internalConnObj=InternalConnectivity(simSpecificInfo);
			%internalConnObj.displayContent();
			
			simSpecificInfoFwd.numPlaces=thisObj.numPlaces;
                        simSpecificInfoFwd.numCellsPerPlace=numCellsPerPlace;
			%thisObj.feedforwardConnObj=copy(FeedForwardConnectivity(simSpecificInfoFwd));
			%thisObj.feedforwardConnObj.displayContent();
			
			thisObj.numCellsL2=Cells.NUM_CELLS_L2;
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%initialize cell state variables
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			initV=NaN(numCellsPerPlace,numPlaces,numSteps);
			initNks=NaN(numCellsPerPlace,numPlaces,numSteps);

			initV(:,:,1)=normrnd(-60,3,numCellsPerPlace,numPlaces);
			initNks(:,:,1)=normrnd(0.1,0,numCellsPerPlace,numPlaces);

			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%group variables into output struct
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			thisObj.numCellsPerPlace=numCellsPerPlace;
			thisObj.numPlaces=numPlaces;
			%{
			thisObj.initV=initV;
			thisObj.initN=initN;
			thisObj.initM=initM;
			thisObj.initH=initH;
			thisObj.initMka=initMka;
			thisObj.initHka=initHka;
			thisObj.initKappaH=initKappaH;
			thisObj.initMnap=initMnap;
			thisObj.initNks=initNks;
			%}
			thisObj.v=initV;
			thisObj.nks=initNks;

			%thisObj.gnapBar=gnapBar;

			%thisObj.gnapBar=gnapBar;
			%thisObj.gnapSigma=gnapSigma;
			thisObj.gksBar=gksBar;
			thisObj.gksSigma=gksSigma;

			%thisObj.injCurrMatrix=injCurrMatrix;

			thisObj.externalInputObj=copy(currInjectorArray);
			thisObj.inhThetaInputArray=copy(inhThetaInputArray);
			thisObj.internalConnObj=copy(internalConnObj);


			thisObj.dt=dt;
			thisObj.numSteps=numSteps;
		end
	end
end


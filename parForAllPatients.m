parpool(10)
%parpool(15)

%parfor ch = 1:96
%parfor ch = 1:192
	%chanNum=ch
	parfor runCount=1:10
	
		 ptSessNum=round(runCount/2)
		stateNum=mod(runCount,2)+1
	%for ptSessNum=13:13
		%for stateNum=1:2
		%for stateNum=1:2
			if(ptSessNum==1)
				ptDir='MG49'
				sessionNum=3
				if(stateNum==1)
					timeLimits=[0 10000];
					stateName='Wake';
				else
					timeLimits=[10000 Inf];
					stateName='Sleep';
				end
			end
			if(ptSessNum==5)
				ptDir='RIHE1'
				sessionNum=13
				if(stateNum==1)
					timeLimits=[7300 19500];
					stateName='Wake';
				else
					continue
				end
			end

			
			if(ptSessNum==2)
				ptDir='MG29'
				sessionNum=2
				if(stateNum==1)
					timeLimits=[-Inf Inf];
					stateName='Wake';
				else
					continue
				end
			end

			
				
			if(ptSessNum==3)
				ptDir='MG29'
				sessionNum=1
				if(stateNum==1)
					timeLimits=[0 9500];
					stateName='Wake';
				else
					timeLimits=[9500 Inf];
					stateName='Sleep';
				end
			end

			
			if(ptSessNum==4)
				ptDir='MG63'
				sessionNum=6
				if(stateNum==1)
					timeLimits=[10000 Inf];
					stateName='Wake';
				else
					timeLimits=[0 10000];
					stateName='Sleep';
				end
			end
			
			
			try
				%doSpectralAnalysesForCh(ptDir,sessionNum,ch,awake1Sleep2,splitStates)
				%doSpectralAnalysesForCh(ptDir,sessionNum,chanNum,timeLimits,stateName)	
				%doSpectralAnalysesForChJustFreqAmpRate(ptDir,sessionNum,chanNum,timeLimits,stateName)	
				%freqPowerDistRatePhasePlot(ptDir,sessionNum,stateName,runCount*2)
				runGetDistantCycleStatsForCell(ptDir,sessionNum,timeLimits,stateName)
			catch ME
				disp(sprintf('error, skipping because:'))
				display(ME.message)
			end
			close all
		end
	%end %pt session loop
%end %ch num loop

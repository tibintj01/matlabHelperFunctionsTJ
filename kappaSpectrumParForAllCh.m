%parpool(19)
parpool(20)

parfor ch = 1:96
%parfor ch = 1:192
	chanNum=ch
	ptDir='MG49'
	sessionNum=3
		%{
		%ptDir='MG63'
		%sessionNum=6
		ptDir='MG29'
		sessionNum=1
		chanNum=ch;
		%}

	if(ch>96)
		sessionNum=2
		chanNum=chanNum-96;	
	end
	for i=1:2
		%awake1Sleep2=i;
		%awake1Sleep2=0;
		%splitStates=0
		if(ch<=96)
			if(i==1)
				timeLimits=[0 10000];
			        stateName='Wake';
				%timeLimits=[0 9500];
			        %stateName='Wake';
				%timeLimits=[10000 Inf];
				%stateName='Wake';
			elseif(i==2)
				timeLimits=[10000 Inf];
				%timeLimits=[9500 Inf];
			        stateName='Sleep';
				%timeLimits=[0 10000];
				%stateName='Sleep-someHighThetaDelta';
			end

		else
			if(i==1)
                                timeLimits=[-Inf Inf];
                                stateName='Wake';
                        elseif(i==2)
                        	continue
			end
		end
		
		try
			%doSpectralAnalysesForCh(ptDir,sessionNum,ch,awake1Sleep2,splitStates)
			doSpectralAnalysesForCh(ptDir,sessionNum,chanNum,timeLimits,stateName)	
		catch ME
			disp(sprintf('error, skipping ch %d because:',chanNum))
			display(ME.message)
		end
	close all
	end
end

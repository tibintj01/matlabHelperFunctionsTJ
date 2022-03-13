for ch=1:96
	try
		checkSeizureSpikeSortingMain(1,ch)
	catch ME
		disp(sprintf('skipping ch %d.....',ch))
		disp(ME.message)
	end
end

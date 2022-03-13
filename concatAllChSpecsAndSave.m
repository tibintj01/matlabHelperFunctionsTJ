
for i=1:96
	if(i==1)
		concatSpecs=getZeroedSgramForCh(i);
	else
		try
			concatSpecs=[concatSpecs getZeroedSgramForCh(i)];
		catch
			disp(sprintf('error, skipping ch %d.....',i))
		end

	end
end

	n=size(concatSpecs,2);
        meanSpec=concatSpecs*ones(n,1)*ones(n,1)'/n;
        concatSpecs=concatSpecs-meanSpec;

A=concatSpecs;

save('MG49-concatenatedSpecs.mat','concatSpecs','A')

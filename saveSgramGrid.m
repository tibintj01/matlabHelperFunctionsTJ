
parfor i=1:22
	chList=[2 68 72 76 80 3 36 40 47 86 6 17 50 56 90 10 23 21 60 94 28 73];
	ch=chList(i)
	try
	getSgramSVDforChNum(ch);
	catch
	disp(sprintf('skipping ch %d',ch))
	end
end

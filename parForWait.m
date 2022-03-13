
%poolobj=gcp
%addAttachedFiles(poolobj,inputPath)
%poolobj

try
	parpool(20)
catch
	display('new worker pool not activated....')
end
tic
parfor i= 1:1000
	A(i)=i;
end
A
toc

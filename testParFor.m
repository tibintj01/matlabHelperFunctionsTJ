disp('serial for.....')
tic
%for i=1:40
%	pause(10)
%end
toc

disp('par for ....')
tic
parfor i=1:100
	pause(60)
end
toc

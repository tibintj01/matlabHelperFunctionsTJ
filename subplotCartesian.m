function [ax] = subplotCartesian(x0,y0,width,height)
	gridSize=100;
	x0Sub=max(1,floor(x0*gridSize));
	x0Sub=min(gridSize,x0Sub);
	y0Sub=max(1,floor(y0*gridSize));
	y0Sub=min(gridSize,y0Sub);

	x1Sub=max(1,floor((x0+width)*gridSize));
	x1Sub=min(gridSize,x1Sub);
	y1Sub=max(1,floor((y0+height)*gridSize));
	y1Sub=min(gridSize,y1Sub);
	
	%subPos1=sub2ind([gridSize gridSize],floor(x0*gridSize), floor(y0*gridSize));
	%subPos2=sub2ind([gridSize gridSize],floor((x0+width)*gridSize), floor((y0+height)*gridSize));
	subPos1=sub2ind([gridSize gridSize],x0Sub,y0Sub);
	subPos2=sub2ind([gridSize gridSize],x1Sub,y1Sub);

	ax=subplot(gridSize,gridSize,[subPos1 subPos2]);

	

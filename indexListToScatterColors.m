function [threeColMatrix]=indexListToScatterColors(idList)

%works for when idList is 1,2,3- extend later!
%scatterCMap= [ 1 0 0;  0 0 1; 0 1 0]; %r,b,g
scatterCMap= [0 0 0;  1 0 0;  0 0 1; 0 1 0; 1 0 1; 1 1 0]; %k,r,b,g,m,y
	
	threeColMatrix=NaN(length(idList),3);

	for i=1:length(idList)
		threeColMatrix(i,:)=scatterCMap(idList(i),:);
	end



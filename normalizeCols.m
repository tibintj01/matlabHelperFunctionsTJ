function [normA]=normalizeCols(A)
	normA=NaN(size(A));
	for i=1:size(A,2)
		currCol=A(:,i);
		normA(:,i)=currCol/norm(currCol);
	end

function [normed] = minNormalizeCols(matr)

	normed=zeros(size(matr));
	for col =1:size(matr,2)
		normed(:,col)=minNormalize(matr(:,col));
	end

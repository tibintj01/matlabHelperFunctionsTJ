function [normed] = rangeNormalizeCols(matr)

	normed=zeros(size(matr));
	for col =1:size(matr,2)
		normed(:,col)=rangeNormalize(matr(:,col));
	end

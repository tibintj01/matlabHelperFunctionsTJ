function [concat]=concatCellArrays(c1,c2)

	c1vert=reshape(c1,[length(c1) 1]);
	c2vert=reshape(c2,[length(c2) 1]);

	concat=[c1vert ; c2vert];

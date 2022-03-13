function [normed]=rangeNormalize(values)
	%normed=scaledata(values,-1,1);	
	%normed=values/(range(values));
	
	%normed=values/(max(values));
	normed=scaledata(values,-1,1);	

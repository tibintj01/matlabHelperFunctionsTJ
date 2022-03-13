function [normed]=rangeNormalize(values)
	%normed=scaledata(values,-1,1);	
	%normed=values/(range(values));
	%normed=-values/(min(values));
	%normed=-values/(min(values));
	normed=scaledata(values,-1,1);	

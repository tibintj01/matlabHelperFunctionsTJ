function [normed]=rangeNormalize(values,prcMin,prcMax)
	%normed=scaledata(values,-1,1);	
	%normed=values/(range(values));
	prcVal1=prctile(values,prcMin);
	prcVal2=prctile(values,prcMax);
	normed=values/(prcVal2-prcVal1);

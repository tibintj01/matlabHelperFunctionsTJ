function [varExplainedFrac]=getVarExplained(sigma,idx)
	varExplainedFrac=sigma(idx)^2/(sum(sigma.^2));

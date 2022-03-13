function [typeChar]=getTypeChar(var)

	if(strcmp(class(var),'char'))
		typeChar='%s';
	elseif(strcmp(class(var),'double'))
		typeChar='%d';
	end

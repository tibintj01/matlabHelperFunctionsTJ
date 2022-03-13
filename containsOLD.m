function [truth]=contains(str,pattern)
	truth=false;
	if(~isempty(findstr(str,pattern)))
		truth=true;
	end

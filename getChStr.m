function [chStr] = getChStr(ch)
	chStr=num2str(ch);
	if(ch<10)
		chStr=['0' chStr];
	end

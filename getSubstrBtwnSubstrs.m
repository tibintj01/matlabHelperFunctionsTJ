function [btwnStr]=getSubstrBtwnSubstrs(str,sub1,sub2)

	k1=strfind(str,sub1)+length(sub1);
	k2=strfind(str,sub2);
	if(length(k2)>0)
		k2=k2(end);
    end

    if(strcmp(sub1,'')) %from start
         btwnStr=str(1:(k2-1));
    elseif(strcmp(sub2,'')) %to end
         btwnStr=str(k1:end);
    else
        btwnStr=str(k1:(k2-1));
    end

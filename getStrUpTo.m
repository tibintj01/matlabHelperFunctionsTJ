function [substr]=getStrUpTo(str,stopChar)

        stopIdxes=findstr(stopChar,str);
        if(isempty(stopIdxes))
            substr=str;
            return
        end
        stopIdx=stopIdxes(1);

        substr=str(1:(stopIdx-1));

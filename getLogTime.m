function [logTimeToEnd] = getLogTime(timeInField,base,noRectify,noLimit)
%properly normalized log transfrom from time to degrees
    %y=(log(x)/log(2)+7)*360/7;
    if(~exist('noLimit') || noLimit==0)
        timeInField(timeInField>1)=NaN;
        timeInField(timeInField<=0)=NaN;
    end
    
   %logTimeToEnd=(log(1-timeInField)/log(2)+7)/7;
    %logTimeToEnd=(log(timeInField)/log(2)+7)/7;
    
    if(~exist('base','var'))
        base=2;
    end
    k=7;
     %k=10;
    logTimeToEnd=(log(timeInField)/log(base)+k)/k;
    %logTimeToEnd=(log(timeInField)/log(base));
       %logTimeToEnd=(log(timeInField)/log(base) + 7);

    %rectify
    if(~exist('noRectify') || noRectify==0)
        logTimeToEnd(logTimeToEnd<0)=0;
    end

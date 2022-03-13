function [logTimeToEnd] = getLogTime(timeInField)
%properly normalized log transfrom from time to degrees
    %y=(log(x)/log(2)+7)*360/7;
    
    onlyInRange=0;
    if(onlyInRange)
        timeInField(timeInField>=1)=NaN;
    end
    timeInField(timeInField<=0)=NaN;
   %logTimeToEnd=(log(1-timeInField)/log(2)+7)/7;
    logTimeToEnd=(log(timeInField)/log(2)+7)/7;
   
    
    if(onlyInRange)
    logTimeToEnd(logTimeToEnd<0)=0;
    end

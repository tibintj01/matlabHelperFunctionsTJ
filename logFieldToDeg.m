function [y] = logFieldFracToDeg(x)
%properly normalized log transfrom from time to degrees
    %y=(log(x)/log(2)+7)*360/7;
   y=(log(x)/log(2)+7)/7;

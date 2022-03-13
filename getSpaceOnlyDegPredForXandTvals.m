function [predictedDegree]=getSpaceOnlyDegPredForXandTvals(xVal,tVal)
%return deg (0, 360) for x (0,1) and t (0,1)
predictedDegree=360*(1-xVal);


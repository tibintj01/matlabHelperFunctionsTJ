function [invertedLogTime] = invertLogTime(logTimeInField,base)

    invertedLogTime=base.^(logTimeInField*7-7);



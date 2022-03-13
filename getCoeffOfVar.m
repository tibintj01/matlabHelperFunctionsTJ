function [cv] = getCoeffOfVar(values)
%returns coefficient of variation

    cv=nanstd(values)/nanmean(values);

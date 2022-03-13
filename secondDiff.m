function [secondDiffVals]=secondDiff(vals)
                secondDiffVals=[0; 0; diff(diff(vals))];

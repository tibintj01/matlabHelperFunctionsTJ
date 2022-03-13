function [firstDiffVals]=firstDiff(vals)
                firstDiffVals=[0; (diff(vals))];

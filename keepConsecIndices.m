function [consecIdxes]=keepConsecIndices(idxes)
 	%keep only consecutive indices
        idxDiffs=[diff(idxes(:));0];
        consecIdxes=idxes(idxDiffs == 1);
end



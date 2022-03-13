function [curveX, curveY] = getBezierGivenConvexity(currConvexityFrac,numPts)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
%if(currConvexityFrac<0)
    if(~exist('numPts','var'))
	numPts=100;
    end

    if(currConvexityFrac<0)
        %[curveX,curveY] = drawBezier([0 -currConvexityFrac 1 1],[1 1 -currConvexityFrac 0],100);
           [curveX,curveY] = drawBezier([0 0 1+currConvexityFrac  1],[0 -currConvexityFrac 1 1],numPts);

    else
        %[curveX,curveY] = drawBezier([0 0 1-currConvexityFrac 1],[1 1-currConvexityFrac 0 0],100);
           [curveX,curveY] = drawBezier([0 currConvexityFrac 1 1],[0 0 1-currConvexityFrac 1],numPts);

    end

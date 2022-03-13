function [rightMinVal,rightMinIdx]=getRightMin(values,binThresh);
        rightValues=values(binThresh:end);
        [rightMinVal,rightMinIdx]=min(rightValues);

	rightMinIdx=rightMinIdx+binThresh-1;

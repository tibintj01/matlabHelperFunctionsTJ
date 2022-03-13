function [leftIdx,leftExtreme,rightIdx,rightExtreme]=getLeftRightExtrema(values)

	[upMaxVal,upMaxIdx]=max(values);
	[downMaxVal,downMaxIdx]=max(-values);


	if(upMaxIdx<downMaxIdx)
		leftExtreme=upMaxVal;
		leftIdx=upMaxIdx;
		rightExtreme=-downMaxVal;
		rightIdx=downMaxIdx;
	else
		leftExtreme=-downMaxVal;
		leftIdx=downMaxIdx;
		rightExtreme=upMaxVal;
		rightIdx=upMaxIdx;
	
	end


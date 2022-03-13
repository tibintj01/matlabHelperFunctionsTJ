function [topIdxes,topValues] = getTopPercentAroundPeak(timeSeries,eventIdx,threshMax)
	%assumes timeSeries contains positive peak at eventIdx, finds first crossing of threshMax fraction 
	%on either side

  maxVal=timeSeries(eventIdx);
  threshVal=maxVal*threshMax;

  firstIdx=max(1,eventIdx-1);
  lastVal=timeSeries(firstIdx);
  while(timeSeries(firstIdx)>threshVal && firstIdx>1)
     lastVal=timeSeries(firstIdx);
     firstIdx=firstIdx-1;
  end

   secondIdx=eventIdx+1;
   lastVal=timeSeries(secondIdx);
   idxLim=length(timeSeries);
   while(secondIdx<idxLim && timeSeries(secondIdx)>threshVal)
      lastVal=timeSeries(secondIdx);
      secondIdx=secondIdx+1;
   end

  topIdxes=firstIdx:secondIdx;
  topValues=timeSeries(firstIdx:secondIdx);

end

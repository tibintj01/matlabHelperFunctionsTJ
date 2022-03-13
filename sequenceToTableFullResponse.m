function [vTrace] = sequenceToTableFullResponse(inputSeq,allInputs,allResponses)
%given sequence (of digits 1-7 and NaN)
%returns subthreshold model response if they were a theta sequence
    sentinelVal=-0.1250;
    startIdx=150;
    startIdx=175;
    
    inputSeq=(inputSeq+0.5)/8;
    inputSeq(inputSeq<0.1)=sentinelVal;
    disp('%%%%%%%%%%%%%%%%%%%%%%%%')
    inputSeq
        disp('%%%%%%%%%%%%%%%%%%%%%%%%')


   %dotProds=allInputs*inputSeq(:);
    %[~,bestMatchIdx]=max(dotProds);
    errVec=sum(abs(allInputs-inputSeq)');
    
    [~,bestMatchIdx]=min(errVec);
    
    vTrace=allResponses(bestMatchIdx,startIdx+1:end);
  %igure
  %plot(vTrace)
  %title(inputSeq)
  %disp('here')
    

    



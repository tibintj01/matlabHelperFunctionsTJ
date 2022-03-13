function [percentOverlap,similarFieldWidths]=getPercentOverlap(currUnitStruct,nextUnitStruct)
%gets quantification of overlap given two fields

similarityToleranceFrac=0.3;
similarFieldWidths=true;

    field1Bounds=[currUnitStruct.fieldPosStart currUnitStruct.fieldPosEnd];
    field1Width=abs(field1Bounds(2)-field1Bounds(1));
    
    field2Bounds=[nextUnitStruct.fieldPosStart nextUnitStruct.fieldPosEnd];
    field2Width=abs(field2Bounds(2)-field2Bounds(1));
    
    avgFieldWidth=(field1Width+field2Width)/2;
    sumNormDiff=abs(field1Width-field2Width)/avgFieldWidth;
    if(sumNormDiff>similarityToleranceFrac)
        similarFieldWidths=false;
    end
     
    allBounds=[field1Bounds field2Bounds];
     
    togetherWidth=abs(min(allBounds(:)) - max(allBounds(:)));
    
    overlapWidth=abs(field2Bounds(1)-field1Bounds(2));
    potentialWidth=field1Width+field2Width;
    
    %percentOverlap=overlapWidth/potentialWidth;
        %percentOverlap=togetherWidth/potentialWidth;
        
        percentOverlap=overlapWidth/togetherWidth;

end


function [dataWithinRect]=getDataWithinRect(data3cols,center,xRad,yRad)

    inXRangeRows=((center(1)+xRad)>data3cols(:,1) & (center(1)-xRad)<data3cols(:,1));
    inYRangeRows=((center(2)+yRad)>data3cols(:,2) & (center(2)-yRad)<data3cols(:,2));
    
    inXandYRangeRows=inXRangeRows & inYRangeRows;
    
    dataWithinRect=data3cols(inXandYRangeRows,:);
    
    %{
	numDataPts=size(data3cols,1);

	dataWithinRect=[];
	for di=1:numDataPts
		currDataPt=data3cols(di,:);

		if(isnan(currDataPt(1)) || isnan(currDataPt(2)))
			continue
		end
		%if x point is outside of center by more than xRad
		if((center(1)+xRad)<currDataPt(1) || (center(1)-xRad)>currDataPt(1))
			continue
		end
		%if y point is outside of center by more than yRad
		if((center(2)+yRad)<currDataPt(2) || (center(2)-yRad)>currDataPt(2))
			continue
		end

		dataWithinRect=[dataWithinRect; currDataPt(:)'];
	end
    %}

	if(isempty(dataWithinRect))
		dataWithinRect=NaN;
	end

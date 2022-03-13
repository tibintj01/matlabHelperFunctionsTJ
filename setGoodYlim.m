function [] = setGoodYlim(values,zeroFloor)

        goodRangeFrac=0.1;
        minY=min(values)-goodRangeFrac*(prctile(values,99)-prctile(values,1));
        maxY=max(values)+goodRangeFrac*(prctile(values,99)-prctile(values,1));
	if(exist('zeroFloor') && zeroFloor==1)
	        minY=max(0, minY)
	end
	ylim([minY maxY])

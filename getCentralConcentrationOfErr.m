function centralConcOfErr=getCentralConcentrationOfErr(fieldTimes,modelErrors)
%returns measure of central concentration of error in time field

maxErrConsidered=100; %deg

    peripheralBound=1/3;
    %peripheralBound=0.2;
    peripheralTimeIdxes=fieldTimes<peripheralBound | fieldTimes>(1-peripheralBound);
    centralTimeIdxes=~peripheralTimeIdxes;
    
    inRangeIdxes=abs(modelErrors)<=maxErrConsidered;
    
    peripheralTimeErrors=modelErrors(peripheralTimeIdxes & inRangeIdxes);
    centralTimeErrors=modelErrors(centralTimeIdxes & inRangeIdxes);


    meanCentralErr=nanmean(centralTimeErrors);
    meanPeripheralErr=nanmean(peripheralTimeErrors);
    %centralConcOfErr=circMeanDeg(centralTimeErrors)-circMeanDeg(peripheralTimeErrors);
        centralConcOfErr=meanCentralErr-meanPeripheralErr;






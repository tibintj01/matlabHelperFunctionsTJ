function [] = saveXCfromCellPropPair(cellPropPath1,cellPropPath2)

	%initialize cross corr settings
	maxLagTime=0.05;
        Fs=2000;
        lags=-maxLagTime*Fs:1:maxLagTime*Fs;
        rsAC=zeros(size(lags));
        fsAC=zeros(size(lags));


	

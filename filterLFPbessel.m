% function [filteredLFP] = filterLFP(lfp, Fs, lowFreq, highFreq,
% filterOrder, zscore)
% Author:        Omar
% Date Created:  5/10/2007
% Date Modified: 11/1/2008
% Notes: Filters LFP, THEN zScores it
% Inputs: filterOrder = order of filter (2 or 20 recommended)
%         filterType = 'butter' or 'equiripple'

function [filteredLFP] = filterLFPbessel(lfp, Fs,lowFreq, filterOrder, zscore)
%This function also used for cell classification code - assumes butterworth filtering
nyqFreq = Fs/2;
%lowFreq = lowFreq/nyqFreq;
%highFreq = highFreq/nyqFreq;

if filterOrder == 6
	%angLowFreq=2*pi*lowFreq;
	%angLowFreq=2*pi*lowFreq/10/10/10/5;
	angLowFreq=2*pi*lowFreq;
    [B, A] = besself(6,angLowFreq,'high' ); % this should be 1, but I have mistakenly been using 2, which means that the default filter I've been using is an order 4 filter; 
    %[B, A] = besself(6,angLowFreq,'low' ); % this should be 1, but I have mistakenly been using 2, which means that the default filter I've been using is an order 4 filter; 
	%this comment is only a concern when passing a third argument of filter type (then uses 2n as order) - Tibin Nov. 2017
    filteredLFP = filtfilt(B, A, double(lfp));

	maxFilteredLFP=max(filteredLFP)

	%%{
	figure
	%maxNumPtsDisp=1e5;
	maxNumPtsDisp=1e4;
	startPt=1
	dispTaxis=linspace(startPt/Fs,maxNumPtsDisp/Fs,maxNumPtsDisp);
	subplot(2,1,1)
	plot(dispTaxis,filteredLFP(startPt:maxNumPtsDisp))
	title(sprintf('filtered bessel high pass %.2f rad/s, %.2f Hz',angLowFreq,angLowFreq/(2*pi)))
	xlim([0.1 Inf])
	subplot(2,1,2)
	plot(dispTaxis,lfp(1:maxNumPtsDisp))
	xlabel('Time (sec)')
	title('raw')
	xlim([0.1 Inf])
	saveas(gcf,sprintf('testBesselHighPass%.3fHz.tif',angLowFreq/(2*pi)))
	%%}
else
    fds %induce runtime error
    %NOT YET IMPLEMENTED
    % Calculate the zpk values using the BUTTER function.
    [z,p,k] = butter(filterOrder/2, [lowFreq highFreq]);
    % To avoid round-off errors, do not use the transfer function.  Instead
    % get the zpk representation and convert it to second-order sections.
    [sos_var,g] = zp2sos(z, p, k);
    hd          = dfilt.df2sos(sos_var, g);
    filteredLFP = filtfilthd(hd, lfp);
end

if zscore == 1
    filteredLFP = zscoreLFP(filteredLFP);
end % end of {if zscore == 1}

clear lfp;

end % end of {function}

function [shiftedVals,bestShift]=xcAlign(vals,trueVals,firstN)
	valsOriginal=vals;
	vals=vals(1:round(firstN));
	trueVals=trueVals(1:round(firstN));
	%vals=vals(30:round(1.5*firstN));
	%trueVals=trueVals(30:(1.5*firstN));

	%convert to col vectors for xcorr!
	vals=vals(:);
	trueVals=trueVals(:);
	

	%trueVals=(trueVals-mean(trueVals))/std(trueVals);
	%vals=(vals-mean(vals))/std(vals);

	%[xc,lag]=xcorr(vals,trueVals,round(firstN))
	%[xc,lag]=xcorr(trueVals,vals,round(firstN))
	%[xc,lag]=xcov(zScore(secondDiff(trueVals)),zScore(secondDiff(vals)),round(firstN))
	%[xc,lag]=xcov(zScore(firstDiff(trueVals)),zScore(firstDiff(vals)),round(firstN))
	[xc,lag]=xcov(zScore(trueVals),zScore(vals),round(firstN));

	[highestCorr,highestCorrIdx]=max(xc);
	bestShift=lag(highestCorrIdx);

	shiftedVals=shiftVals(valsOriginal,bestShift);
	%close all; plot(vals)
	%hold on; plot(trueVals,'k')
	%plot(shiftedVals,'r')
	%saveas(gcf,'testXCalign.tif')

	%fds	
	

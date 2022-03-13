cellPropPaths=getRegexFilePaths('/nfs/turbo/lsa-ojahmed/classifiedMG49-sleepWake','*cell_properties*mat');
close all
rsCount=0;
fsCount=0;
figure(1)
rsHa = tight_subplot(9,9,[.01 .03],[.1 .01],[.01 .01]);

figure(2)
fsHa = tight_subplot(9,9,[.01 .03],[.1 .01],[.01 .01]);

for i=1:length(cellPropPaths)
	[rsCount,fsCount]=plotSpikeVsNonSpikeLFPDist(cellPropPaths{i},rsCount,fsCount,rsHa,fsHa);

end

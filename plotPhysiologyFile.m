function [] = plotPhysiologyFile(filePath)

	[data,sampPeriod,header]=abfload(filePath);
	microSecToSec=1/1000000;	

	singleSweep=squeeze(data(:,2,10));
	timeAxis=(1:length(singleSweep))*sampPeriod*microSecToSec;


	figure
	plot(timeAxis,singleSweep,'k')
	xlim([0 2])
	ylabel('Transmembrane potential (mV)')
	xlabel('Time (sec)')

	plotVertLine(1,'b')
	plotVertLine(1.6,'r')
	legend('Targeted cell activity', 'Light onset','Light offset','Location','Best')
	set(gca,'FontSize',16)


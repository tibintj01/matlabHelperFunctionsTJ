

%filePaths=getRegexFilePaths('/home/tibintj/turboHome/spikeDynamicsAnalysisTibin/compiledDir_02-Nov-2017_04-02-52_getCenterOfMassBinsPar','center*mat')
filePaths=getRegexFilePaths('/home/tibintj/turboHome/spikeDynamicsAnalysisTibin/compiledDir_02-Nov-2017_14-35-31_getCenterOfMassBinsPar','center*mat')
close all
vals=zeros(size(filePaths));

for i=1:length(filePaths)
	i
	data=matfile(filePaths{i});
	%vals(i)=-(data.troughCMBin - data.peakCMBin);
	plot(data.troughCMBin, data.peakCMBin,'ko','MarkerSize',1)
	hold on
end
xlabel('troughCMbin')
ylabel('PeakCMbin')
%histogram(vals)

saveas(gcf,'troughCMbinVsPeakCMbinAvgWaveformsB4Sz.tif')



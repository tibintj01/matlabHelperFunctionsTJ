function [] = saveAsTif(fileName) 
	%saveas(gcf,sprintf('%s.tif',fileName));
	print('-r600',gcf,fileName,'-dtiff')

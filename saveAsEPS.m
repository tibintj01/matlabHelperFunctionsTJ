function [] = saveAsEPS(fileName,dpiVal,onlySVG) 
	saveas(gcf,sprintf('%s.tif',fileName));

	if(exist('dpiVal') && dpiVal>0)
		dpiStr=sprintf('-r%d',dpiVal);
	else
		dpiStr='-r300';
	end

	%figure('Renderer','Painters')
	if(exist('onlySVG') && onlySVG==1)
		saveFileFormatStr='-dsvg';
	else
		saveFileFormatStr='-depsc';
		%saveFileFormatStr='-depsc2';
	end

	%print(dpiStr,gcf,fileName,'-depsc','-tiff')
	print(dpiStr,gcf,fileName,saveFileFormatStr,'-tiff')

	

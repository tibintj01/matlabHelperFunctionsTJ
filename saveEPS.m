function [] = saveAsEPS(filePath,dpiVal) 


if(~contains(filePath,'.eps'))
    filePath=[filePath '.eps'];
        
end
	if(exist('dpiVal') && dpiVal>0)
		dpiStr=sprintf('-r%d',dpiVal);
	else
		dpiStr='-r300';
    end


 saveFileFormatStr='-depsc';



	print(dpiStr,gcf,filePath,saveFileFormatStr)

	

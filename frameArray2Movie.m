function [] = frameArray2Movie(savePath,F)

	writerObj=VideoWriter(savePath);
	writerObj.FrameRate=10;

	open(writerObj)

	for i=1:length(F)
		frame=F(i);
		writeVideo(writerObj,frame);
	end

	close(writerObj);


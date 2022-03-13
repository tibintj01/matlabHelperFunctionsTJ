function [] = matrix2MovieTibin(fileName,matrix3D,frameRate,colormap)

	%add color dimension
	if(length(size(matrix3D))==3)
			
		matrix4D=NaN(size(matrix3D,1),size(matrix3D,2),1,size(matrix3D,3));
		for i=1:size(matrix3D,3)
			matrix4D(:,:,1,i)=matrix3D(:,:,i);
		end
		matrix3D=matrix4D;
		matrix3D=discretize(matrix3D,length(colormap));
	end

	mov=immovie(matrix3D,colormap);
	%implay(mov)

	disp('saving movie.....')
	v=VideoWriter(fileName)
	v.FrameRate=frameRate;
	open(v)
	%for k=1:size(matrix3D,4)
	%	frame=uint8(matrix3D(:,:,1,k));
	%	writeVideo(v,frame)
	%end
	writeVideo(v,mov)
	close(v)



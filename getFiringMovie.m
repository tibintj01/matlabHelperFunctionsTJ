function [] = getFiringMovies(dirName)
	cellPropFilePaths=getFilePathsRegex(dirName,'*.mat');
	tic
	maxTime=300;

	rateWind=5;
	rs3D=zeros(10,10,1,100000);
	fs3D=zeros(10,10,1,100000);
	disp('populating 3D matrix of firing rates by cell type.....')
	%save firing rate distributions by cell type	
	maxLengthSoFar=0;
	for cellIdx=1:length(cellPropFilePaths)
		cellProp=load(cellPropFilePaths{cellIdx});
		
		firingRate=getFiringRateOverTime(cellProp.spikeTimes,rateWind);

		[row,col,cellIdx]=getSpatialRowColFromFile(cellPropFilePaths{cellIdx});
		if(cellProp.isInterneuronCell==0)
			rs3D(row,col,1,1:length(firingRate))=firingRate;
		elseif(cellProp.isInterneuronCell==1)
			fs3D(row,col,1,1:length(firingRate))=firingRate;
		end

		if(length(firingRate)>maxLengthSoFar)
			maxLengthSoFar=length(firingRate);
		end
	end	
	toc

	rs3D=scaledata(rs3D,1,255);
	fs3D=scaledata(fs3D,1,255);

	stopFrame=min(maxLengthSoFar,round(maxTime/rateWind))
	rs3D=rs3D(:,:,1,1:stopFrame);
	fs3D=fs3D(:,:,1,1:stopFrame);
	%rs3D=rs3D(:,:,1,1:maxLengthSoFar);
	%fs3D=fs3D(:,:,1,1:maxLengthSoFar);


	rs3D=repelem(rs3D,10,10);
	fs3D=repelem(fs3D,10,10);


	%fds

	%rs3D=uint16(rs3D);
	%fs3D=uint16(fs3D);
	%rs3D=rs3D(:,:,:,1:100);
	%fs3D=fs3D(:,:,:,1:100);

	frameRate=1/rateWind;
	matrix2MovieTibin(fullfile(dirName,sprintf('RS-SpatialFiring-%.2e-secWind.avi',rateWind)),rs3D,frameRate,jet)
	matrix2MovieTibin(fullfile(dirName,sprintf('FS-SpatialFiring-%.2e-secWind.avi',rateWind)),fs3D,frameRate,jet)



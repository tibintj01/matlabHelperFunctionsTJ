function []= plotFeatures2D(featureTable,featureErrTable,fN1,fN2,groupName,isInterneuron,saveDir,Fs)
	
ephysFigureColorsRatRSC
plotSEMs=0;
%figure;plot3(featureTable{:,butter300T15,featureTable.nexTPwidth,featureTable.nexP50,'ko','MarkerSize',2)
	if(~exist('isInterneuron','var'))
		figure;plot(featureTable{:,fN1},featureTable{:,fN2},'ko','MarkerSize',32)
	else
		%red and blue
		colorVec=repmat(RSColor,height(featureTable),1)
		for i=1:length(colorVec)
			%if(isnan(isInterneuron(i)) || isInterneuron(i)==2)
			if(isInterneuron(i)==2)
				%colorVec(i,:)=[0 0 0];%FS1ColorLight;
				colorVec(i,:)=FS2Color;
			elseif (isInterneuron(i))
				colorVec(i,:)=FS1Color;
			end
		end
		%figure;scatter3(featureTable{:,fN1},featureTable{:,fN2},featureTable{:,fN3},2,colorVec)
		if(exist('Fs','var')) %only if plotting temporal feature...
			figure;
			%scatter(featureTable{:,fN1}/Fs*1000,featureTable{:,fN2}/Fs*1000,32,colorVec)
		
			if(contains(fN1,'atio') ||  contains(fN1,'Amp'))
				scatter(featureTable{:,fN1},featureTable{:,fN2}/Fs*1000,1,colorVec)
				xlabel(sprintf('%s',fN1))
				ylabel(sprintf('%s (ms)',fN2))
			elseif(contains(fN2,'atio') ||  contains(fN2,'Amp'))
				scatter(featureTable{:,fN1}/Fs*1000,featureTable{:,fN2},1,colorVec)
				ylabel(sprintf('%s',fN1))
				xlabel(sprintf('%s (ms)',fN2))
			else
				scatter(featureTable{:,fN1}/Fs*1000,featureTable{:,fN2}/Fs*1000,1,colorVec)
				xlabel(sprintf('%s (ms)',fN1))
				ylabel(sprintf('%s (ms)',fN2))
			end
			%hold on
			
			if(plotSEMs)
				xErr=featureErrTable{isInterneuron==0,fN1}/Fs*1000;
				yErr=featureErrTable{isInterneuron==0,fN2}/Fs*1000;

				%errorbar(featureTable{isInterneuron==0,fN1}/Fs*1000,featureTable{:,fN2}/Fs*1000,yErr,yErr,xErr,xErr,'o','MarkerSize',1,')			
				errorbar(featureTable{isInterneuron==0,fN1}/Fs*1000,featureTable{isInterneuron==0,fN2}/Fs*1000,yErr,yErr,xErr,xErr,'.r','CapSize',4)
				hold on			
				
				xErr=featureErrTable{isInterneuron==1,fN1}/Fs*1000;
				yErr=featureErrTable{isInterneuron==1,fN2}/Fs*1000;
				errorbar(featureTable{isInterneuron==1,fN1}/Fs*1000,featureTable{isInterneuron==1,fN2}/Fs*1000,yErr,yErr,xErr,xErr,'.b','CapSize',4)			
				
				xErr=featureErrTable{isInterneuron==2,fN1}/Fs*1000;
				yErr=featureErrTable{isInterneuron==2,fN2}/Fs*1000;
				errorbar(featureTable{isInterneuron==2,fN1}/Fs*1000,featureTable{isInterneuron==2,fN2}/Fs*1000,yErr,yErr,xErr,xErr,'.k','CapSize',4)			
		
			end
		else
			figure;scatter(featureTable{:,fN1},featureTable{:,fN2},32,colorVec)
			xlabel(fN1)
			ylabel(fN2)
		end
	end
	title(sprintf('%s, %d total',removeUnderscores(groupName),length(isInterneuron)))	
	if(~exist('saveDir','var'))
		saveas(gcf,sprintf('%s-2D-%svs%s.tif',groupName,fN1,fN2))
	else
		saveas(gcf,fullfile(saveDir,sprintf('%s-2D-%svs%s.tif',groupName,fN1,fN2)))
	end

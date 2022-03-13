function []= plotFeatures2D(featureTable,fN1,fN2,groupName,isInterneuron,saveDir,Fs,ax)
	
ephysFigureColors

%if(exist('ax','var'))
%	axes(ax)
%end

%figure;plot3(featureTable{:,butter300T15,featureTable.nexTPwidth,featureTable.nexP50,'ko','MarkerSize',2)
	if(~exist('isInterneuron','var'))
		plot(featureTable{:,fN1},featureTable{:,fN2},'ko','MarkerSize',32)
	else
		%red and blue
		colorVec=repmat(RSColor,height(featureTable),1)
		for i=1:length(colorVec)
			if(isnan(isInterneuron(i)) || isInterneuron(i)==2)
				colorVec(i,:)=[0 0 0];
			elseif (isInterneuron(i))
				colorVec(i,:)=FS1Color;
			end
		end
		%figure;scatter3(featureTable{:,fN1},featureTable{:,fN2},featureTable{:,fN3},2,colorVec)
		if(exist('Fs','var')) %only if plotting temporal feature...
			%scatter(featureTable{:,fN1}/Fs*1000,featureTable{:,fN2}/Fs*1000,32,colorVec)
			scatter(featureTable{:,fN1}/Fs*1000,featureTable{:,fN2}/Fs*1000,2,colorVec)
			%xlabel(sprintf('%s (ms)',fN1))
			%ylabel(sprintf('%s (ms)',fN2))
			xlabel('15% Trough Width (ms)')
			ylabel(sprintf('Trough to Peak Width (ms)'))
		else
			%scatter(featureTable{:,fN1},featureTable{:,fN2},32,colorVec)
			scatter(featureTable{:,fN1},featureTable{:,fN2},2,colorVec)
			xlabel(fN1)
			ylabel(fN2)
		end
	end
	title(removeUnderscores(groupName))	
			pbaspect([3 4 1])	
			print('-r600',gcf,'ClassificationWaveformFeatureSpace','-depsc','-tiff')	
%	if(~exist('ax','var'))
		if(~exist('saveDir','var'))
			pbaspect([3 4 1])	
			print('-r600',gcf,'ClassificationWaveformFeatureSpace','-depsc','-tiff')	
			saveas(gcf,sprintf('%s-2D-%svs%s.tif',groupName,fN1,fN2))
		else
			saveas(gcf,fullfile(saveDir,sprintf('%s-2D-%svs%s.tif',groupName,fN1,fN2)))
			
		end
%	end

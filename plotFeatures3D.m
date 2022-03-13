function []= plotFeatures3D(featureTable,fN1,fN2,fN3,groupName,isInterneuron)
	

%figure;plot3(featureTable{:,butter300T15,featureTable.nexTPwidth,featureTable.nexP50,'ko','MarkerSize',2)
	if(~exist('isInterneuron','var'))
		figure;plot3(featureTable{:,fN1},featureTable{:,fN2},featureTable{:,fN3},'ko','MarkerSize',2)
	else
		%red and blue
		colorVec=repmat([1 0 0],height(featureTable),1)
		for i=1:length(colorVec)
			if(isnan(isInterneuron(i)))
				colorVec(i,:)=[0 0 0];
			elseif (isInterneuron(i))
				colorVec(i,:)=[0 0 1];
			end
		end
		figure;scatter3(featureTable{:,fN1},featureTable{:,fN2},featureTable{:,fN3},2,colorVec)
	end
	xlabel(fN1)
	ylabel(fN2)
	zlabel(fN3)
	
	saveas(gcf,sprintf('%s-3D-%svs%svs%s.tif',groupName,fN1,fN2,fN3))


function [inputs] = interactiveSVDspace(U,V,k,timeAxis)
%makes interactive plot with first k left singular vectors

	fHandle=@getUVecComb;
	initComb=zeros(k,1);
	initComb(1)=1;

	maxY=max(max(U(:,1:k)))*1.5;
	minY=min(min(U(:,1:k)))*1.5;

	vecU=U(:,1:k);
	vecU=vecU(:);

	minY=-0.5;
	maxY=0.5;

	%inputs={{fHandle,vecU,initComb},{[min(timeAxis) max(timeAxis)],[minY maxY]}};
	inputs={{fHandle,timeAxis,initComb},{[min(timeAxis) max(timeAxis)],[minY maxY]}};

	for i=1:k
		inputs{2+i}={i,sprintf('Comp. %d coeff',i),min(V(:,i)),max(V(:,i))};
	end
	%inputs
	manipulate(inputs{:})


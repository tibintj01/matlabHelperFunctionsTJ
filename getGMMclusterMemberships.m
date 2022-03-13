function [isInterneuron,p,M]=getGMMclusterMemberships(featureMatrix,numClusters,grayPVal,ax)

	ephysFigureColors

	plotOriginalClusters=0;

	options = statset('Display','final');
	gm=fitgmdist(featureMatrix,numClusters,'Options',options,'Replicates',1000)
	[idx,nlogl,P,logpdf,M]=cluster(gm,featureMatrix);

	
	P=posterior(gm,featureMatrix);
	plotProb=1;
	swapClusters=0;
	
	%posterior probablity of being in 1st cluster given that this data was observed
	pInRSgauss=P(:,1);
	pInFSgauss=P(:,2);

	if(plotProb)
	%adapted from https://www.mathworks.com/help/stats/cluster-data-from-mixture-of-gaussian-distributions.html
		cluster1 = (idx == 1); % |1| for cluster 1 membership
		cluster2 = (idx == 2); % |2| for cluster 2 membership

		%determine if cluster1 is RS or FS (which has the greater mean widths) - cluster1 should be RS; perhaps not
		%if(nanmean(featureMatrix(cluster1,2))<nanmean(featureMatrix(cluster2,2))) - is 1 or 2 RS vs FS, where is this determined?
		if(nanmedian(featureMatrix(cluster1,2))<nanmedian(featureMatrix(cluster2,2)) || nanmedian(featureMatrix(cluster1,1))<nanmedian(featureMatrix(cluster2,1)))
             		%removing swapping makes right number of FS?
			%temp=cluster2;
                    %cluster2=cluster1;
                    %cluster1=temp;
			swapClusters=1;
        	end
	

		if(swapClusters)
			temp=pInRSgauss;
			pInRSgauss=pInFSgauss;
			pInFSgauss=temp;
			%probIsRS=(1-p);
			M(:,[1 2])= M(:,[2 1]); 
		end

		if(plotOriginalClusters)
			figure;
			%scatter(featureMatrix(cluster1,1),featureMatrix(cluster1,2),15,P(cluster1,1),'+')
			scatter(featureMatrix(cluster1,1),featureMatrix(cluster1,2),2,P(cluster1,1),'o')
			hold on
			%scatter(featureMatrix(cluster2,1),featureMatrix(cluster2,2),15,P(cluster2,1),'o')
			scatter(featureMatrix(cluster2,1),featureMatrix(cluster2,2),2,P(cluster2,1),'o')
			%set(gca,'Color','k')
			hold off
			%clrmap = jet(80);
			%clrmap = copper(80);
			%colormap(clrmap(9:72,:))
			colormap(bkgColorMap)
		%xlabel('15% Trough width')
		%ylabel('Trough to peak width')
			ylabel(colorbar,'RS Cluster Posterior Prob.')
			title('RS Cluster Posterior Probabilites')
			saveas(gcf,'originalClustering.tif')		
		end

		if(exist('ax','var'))
			%axes(ax)
			figure
		else
			figure
		end
	        
		scatter(M(cluster1,1),M(cluster1,2),2,P(cluster1,1),'o')
                hold on
                scatter(M(cluster2,1),M(cluster2,2),2,P(cluster2,1),'o')
                hold off	
		%clrmap = copper(80);
		%colormap(clrmap(9:72,:))
		colormap(bkgColorMap)
		xlabel('M distance from RS cluster')
		ylabel('M distance from FS cluster')
		ylabel(colorbar,'RS Cluster Posterior Prob.','FontSize',10)
		title('Mahalanobis distances')
		pbaspect([3 4 1])
		print('-r600',gcf,'ClassifiedUnit_Mdistances','-depsc','-tiff')
		%saveas(gcf,'classifyTest.tif')
	end

	

		
	pThresh=1-grayPVal;

	%isInterneuron=zeros(size(featureMatrix,1),1);
	isInterneuron=2*ones(size(featureMatrix,1),1);
	%isInterneuron(pInRSgauss<grayPVal)=1;
	isInterneuron(pInFSgauss>pThresh)=1;
	isInterneuron(pInRSgauss>pThresh)=0;

	p=pInRSgauss;	

	%cutoff of posterior for intermediate with RS? 95
	%isInterneuron(pInRSgauss>grayPVal & pInRSgauss< (1-grayPVal))=2; %2 signifies gray
	
	

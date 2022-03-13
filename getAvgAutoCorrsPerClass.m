function [] = getAvgAutoCorrsPerClass(cellPropDir)
	cellPropFilePaths=getFilePathsRegex(cellPropDir,'*.mat');
        tic
     
	rsCount=0;
	fsCount=0; 

	%these should match getCrossCorr.m.....
	maxLagTime=0.05;
	Fs=2000;
	lags=-maxLagTime*Fs:1:maxLagTime*Fs;
	rsAC=zeros(size(lags));
	fsAC=zeros(size(lags));

	figure
	numFSSpikes=0;
	numRSSpikes=0;
        for cellIdx=1:length(cellPropFilePaths)
                cellProp=load(cellPropFilePaths{cellIdx});

		%if(rsCount>100 && cellProp.isInterneuronCell==0)
		%	continue
		%end
      		[xc,normxc]=getCrossCorr(cellProp.spikeTimes(cellProp.goodSpikes),cellProp.spikeTimes(cellProp.goodSpikes),'test',0);
		numSpikes=length(cellProp.spikeTimes(cellProp.goodSpikes));

		if(numSpikes>0 && isempty(find(isnan(normxc))))

			if(cellProp.isInterneuronCell==0)
				rsCount=rsCount+1;
				rsAC=rsAC+normxc;		
				numRSSpikes=numRSSpikes+length(cellProp.spikeTimes(cellProp.goodSpikes));
			elseif(cellProp.isInterneuronCell==1)
				fsCount=fsCount+1;
				fsAC=fsAC+normxc;		  
				numFSSpikes=numFSSpikes+length(cellProp.spikeTimes(cellProp.goodSpikes));
		
			end
		end
	hold off
	plot(lags/Fs*1000,rsAC/rsCount,'r')
	hold on
	plot(lags/Fs*1000,fsAC/fsCount,'b')
	legend('RS','FS')
	title('Spike train autocorrelation')

                xlabel('Spike train lag (msec)')
                ylabel('Overlapping spike count (normalized)')
	
	drawnow
	end
	

	saveas(gcf,'avgAutoCorrelationsByClass.tif')	

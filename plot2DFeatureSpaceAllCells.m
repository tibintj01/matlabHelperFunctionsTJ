function [] = plot2DFeatureSpaceAllCells(featureIdx1,featureIdx2,titleSuffix)
	if(isstr(featureIdx1))
		featureIdx1=str2num(featureIdx1);
	end
	if(isstr(featureIdx2))
		featureIdx2=str2num(featureIdx2);
	end
	%try

		%cellPropsDir='/home/tibintj/seizureCellWaveformProperties';
		%cellPropsDir='/home/tibintj/compiledDir_01-Nov-2017_18-50-18_getCellWaveformPropsOmar_Tibin/seizureCellWaveformProperties';
		%cellPropsDir='/home/tibintj/compiledDir_05-Nov-2017_01-37-49_getCellWaveformPropsOmar_Tibin/seizureCellWaveformProperties-05-Nov-2017_01-36-44'
		%raw waveform properties
		%cellPropsDir='/home/tibintj/pbsTest/compiledDir_17-Nov-2017_13-07-44_getCellWaveformPropsOmar_Tibin/seizureCellWaveformProperties-05-Nov-2017_01-36-44';
		%cellPropsDir='/home/tibintj/compiledDir_19-Nov-2017_02-27-27_getCellWaveformPropsOmar_Tibin/seizureCellWaveformProperties-05-Nov-2017_01-36-44';
		%cellPropsDir='/home/tibintj/compiledDir_20-Nov-2017_19-50-56_getCellWaveformPropsOmar_Tibin/seizureCellWaveformPropertiesRaw-20-Nov-2017_19-07-32';
		%cellPropsDir='/home/tibintj/compiledDir_21-Nov-2017_01-52-38_getCellWaveformPropsOmar_Tibin_sepSeizureNEXFilt/seizureCellWaveformPropertiesNEX-21-Nov-2017_01-52-34';
		%cellPropsDir='/home/tibintj/compiledDir_21-Nov-2017_03-45-10_getCellWaveformPropsOmar_Tibin_sepSeizureButter300Filt/seizureCellWaveformPropertiesNEX-21-Nov-2017_01-52-34/';
		cellPropsDir='/home/tibintj/turboHome/spikeDynamicsAnalysisTibin/compiledDir_22-Nov-2017_01-17-18_getCellWaveformPropsOmar_Tibin_Raw-correctedCellNo/seizureCellWaveformProperties-22-Nov-2017_00-32-58';



		cellPropPaths=getRegexFilePaths(cellPropsDir,'*cell_properties*.mat');

		featureIdxOffset=4;
		figIdxes=[];

		close all
		numCells=length(cellPropPaths);
		for cellIdx=1:numCells
			if(mod(cellIdx,50)==0 || cellIdx==numCells)
				disp(sprintf('Getting features for cell %d....',cellIdx))
			end
			cellProps=load(cellPropPaths{cellIdx});
		
			[feature1Name,feature1array]=getNthField(cellProps,featureIdx1+featureIdxOffset);
			[feature2Name,feature2array]=getNthField(cellProps,featureIdx2+featureIdxOffset);
		

			%szStartTimes=getSeizureStartTimes(cellPropPaths{cellIdx});
			%szStartTime=szStartTimes(1);

			%szColor=getSeizureColor(cellPropPaths{cellIdx});

			%if(length(feature1array)>1)
			%	feature1arrayB4sz=feature1array(cellProps.spikeTimes<szStartTime);	
			%else
			%	feature1arrayB4sz=feature1array;
			%end
		
			%if(length(feature2array)>1)
			%	feature2arrayB4sz=feature2array(cellProps.spikeTimes<szStartTime);	
			%else
			%	feature2arrayB4sz=feature2array;
			%end

			%feature1Val=nanmean(feature1arrayB4sz);
			%feature2Val=nanmean(feature2arrayB4sz);
			feature1Val=nanmean(feature1array);
			feature2Val=nanmean(feature2array);

			figIdx=str2num(sprintf('%d%d',getSessionIdx(cellPropPaths{cellIdx}),getSzIdx(cellPropPaths{cellIdx})));	
			%if index not in list of figure indices, add it to list
			if(isempty(figIdxes) || min(abs(figIdxes-figIdx))>0.01)
				figIdxes=[figIdxes figIdx];
				figureTitles{length(figIdxes)}=sprintf('%s-Seizure%d',getSessionName(cellPropPaths{cellIdx}),getSzIdx(cellPropPaths{cellIdx}));
			end

			figure(figIdx)		
			if(getSessionIdx(cellPropPaths{cellIdx})==1)
			%	fds
			end
			%plot(feature1Val,feature2Val,[szColor 'o'],'MarkerSize',1)
			plot(feature1Val,feature2Val,['k' 'o'],'MarkerSize',1)
			hold on
		end
		
		%if(featureIdx1>18 && featureIdx2>18)
		%	daspect([1 1 1])
		%	hold on
			%plot([min(feature1Val)*1.2 max(feature1Val)],[0 0],'k--')
		%	plot([-1 1],[0 0],'k--')
		%	xlim([-0.9 0.5])
		%	ylim([-0.8 0.3])
			%xlim([min(feature1Val)*1.2 max(feature1Val)])
			%plot([0 0],[min(feature2Val)*1.2 max(feature2Val)],'k--')	
		%	plot([0 0],[-1 1],'k--')
			%ylim([min(feature2Val)*1.2 max(feature2Val)])
		%end

		%[colors,names]=getSeizureColorNames();
		%addCustomLegendTibin(colors,names)

		%title('All seizure session cells before seizure')
		outDirName=['/nfs/turbo/lsa-ojahmed/feature2DplotsOutputTibinParEachSeizure-' titleSuffix] ;
		if(~isdir(outDirName))
			mkdir(outDirName)
		end

		for i=1:length(figIdxes)
			figure(figIdxes(i))
			title([figureTitles{i} '-' titleSuffix])
			xlabel(feature1Name)
			ylabel(feature2Name)
			saveas(gcf,sprintf('%s/seizureAllCells-%s-Tibin-%s-vs-%s.tif',outDirName,figureTitles{i},feature1Name,feature2Name))
		end

	%catch
	%	disp(sprintf('Error for %d vs %d! Did not run.....',featureIdx1,featureIdx2))
	%end

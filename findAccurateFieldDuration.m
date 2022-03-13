close all; clear all; clc
dataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles';
showPlots=0;
showPlots=1;

filePaths=getFilePathsRegex(dataDir,'*mat');

for i=1:length(filePaths)
      i/length(filePaths)
    currDataFilePath=filePaths{i};
    
    currFileName=getFileNameFromPath(currDataFilePath);
    fileBaseName=currFileName(1:(end-4));
    
      currData=load(currDataFilePath);
      
      for di=1:2
          if(di==2)
            currFieldDirection='leftward';  
          else
              currFieldDirection='rightward';  
          end
          
          if(~isfield(currData.directionSpecificStats,currFieldDirection))
              continue
          end
          
          numFields=size(currData.directionSpecificStats.(currFieldDirection).fieldWidthsPerLap,1);
          
          for fi=1:numFields
              currFieldDurationPerLap=currData.directionSpecificStats.(currFieldDirection).fieldDurationsPerLap(fi,:);
              currFieldWidthPerLap=abs(currData.directionSpecificStats.(currFieldDirection).fieldWidthsPerLap(fi,:));
              
              lapInterferenceBool=currData.hasLapFieldInterference;
              if(di>size(lapInterferenceBool,1) || fi > size(lapInterferenceBool,2))
                  continue
              end
              currLapInterferenceBool=logical(squeeze(lapInterferenceBool(di,fi,:)));
              
              currFieldDurationPerLap(currLapInterferenceBool)=NaN;
              currFieldWidthPerLap(currLapInterferenceBool)=NaN;
              
              currFieldAvgSpeedPerLap=currFieldWidthPerLap./currFieldDurationPerLap;
              nBins=20;
              widthEdges=0:0.01:1;
              durationEdges=0:0.1:10;
              speedEdges=0:0.01:1;

              if(sum(isnan(currFieldWidthPerLap))==length(currFieldWidthPerLap))
                  continue
              end

              if(sum(isnan(currFieldAvgSpeedPerLap))==length(currFieldAvgSpeedPerLap))
                  continue
              end

              figure;
              subplot(2,4,1); histogram(currFieldWidthPerLap,widthEdges); xlabel('field width per lap (m)'); ylabel('lap count')

              subplot(2,4,2); histogram(currFieldDurationPerLap,durationEdges);xlabel('field duration per lap (sec)'); ylabel('lap count')
              subplot(2,4,3); histogram(currFieldAvgSpeedPerLap,speedEdges);xlabel('avg field speed per lap (m/sec)'); ylabel('lap count')
              subplot(2,4,4); plot(currFieldWidthPerLap,currFieldDurationPerLap,'k.'); xlabel('field width per lap (m)');ylabel('field duration per lap (sec)');
              xlim([0 1])
              ylim([0 5])

              if(~iscell(currData.imagePaths))
                  if(~contains(currData.imagePaths,currFieldDirection))
                      continue
                  end
                subplot(2,4,[ 5 6 7 8]); imshow(currData.imagePaths)
              else
                  if(di==1)
                      ii=2;
                  else
                      ii=1;
                  end
                   if(~contains(currData.imagePaths{ii},currFieldDirection))
                      continue
                  end
                  subplot(2,3,[4 5 6]); imshow(currData.imagePaths{ii})
              end

              uberTitle(removeUnderscores(sprintf('%s %s, field %d',fileBaseName,currFieldDirection,fi)))
              maxFig
              setFigFontTo(18)

              saveas(gcf,sprintf('%s_FieldWidthDurationStats_%s_field%d.tif',fileBaseName,currFieldDirection,fi))
              close all
          end%field loop
      end
end

  
dataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles';

sessionNames={'Achilles_10252013','Achilles_11012013','Buddy_06272013','Cicero_09012014','Cicero_09102014','Cicero_09172014','Gatsby_08282013'};

numChannels=128;

for si=1:length(sessionNames)
    currSessName=sessionNames{si};
    
    if(strcmp(currSessName,'Cicero_09012014') || strcmp(currSessName,'Cicero_09172014'))
        continue
    end
    
    filePaths=getFilePathsRegex(dataDir,sprintf('%s*mat',currSessName));
    numUnits=length(filePaths);
    allNormChLocking=NaN(numChannels,numUnits);
    for fi=1:numUnits
        currUnitData=load(filePaths{fi});
        allNormChLocking(:,fi)=currUnitData.chKappasZscore(:);
        
    end
    meanNormChLocking=nanmean(allNormChLocking,2);
    semNormChLocking=getSEMacrossRows(allNormChLocking');
    figure;
    pH=plot(allNormChLocking,'k.');
    %alpha(0.1)
    hold on
    shadedErrorBar(1:numChannels,meanNormChLocking,semNormChLocking)
    %axis tight
    ylim([-2 2])
    xlim([1 numChannels])
    xticks(1:128)
    grid on
    %xticklabels({1:10:120})
    xlabel('Ch number')
    ylabel('Phase locking')
    legend(sprintf('n=%d units',numUnits))
    [~,maxLockingCh]=max(meanNormChLocking);
    title(removeUnderscores(sprintf('%s, max locking channel: %d',currSessName, maxLockingCh)))
  
     setFigFontTo(18)
   a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',7)
 
    %title(removeUnderscores(currSessName))
     maxFig
      
    saveas(gcf,sprintf('%s_PhaseLockingAcrossChannelsAboveMinSpeed.tif',currSessName))
    close all
end

close all; clear all; clc
%recalculate=0;

if(~exist('durationStatsOct14.mat','file'))
    processedDataDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData';
    unitDataDir=fullfile(processedDataDir,'unitSpikeInfoStructs');

    %savePhaseVsThetaCyclesElapsedDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3Images/phaseVsTimeInFieldInCycles';
    %touchDir(savePhaseVsThetaCyclesElapsedDir)

    filePaths=getFilePathsRegex(unitDataDir,'*mat');


    startFile=1;

    allCycleTimeBounds=[];

    allNormFieldDurationsPerLap=[];
    allAbsoluteFieldDurationsPerLap=[];
    dirStr={'rightward','leftward'};

    for fi=startFile:length(filePaths)

        fi
        currFilePath=filePaths{fi};
        data=load(currFilePath);

        rightCycleTimeBounds=data.cycleTimeFieldBoundPerDirPerField.rightward;
        leftCycleTimeBounds=data.cycleTimeFieldBoundPerDirPerField.leftward;

        allCycleTimeBounds=[allCycleTimeBounds;rightCycleTimeBounds;leftCycleTimeBounds];

        for di=1:2

            if(~isfield(data,'thetaBasedTimeVsPhaseInfo'))
                continue
            end

            if(~isfield(data.thetaBasedTimeVsPhaseInfo,dirStr{di}))
                continue
            end

            currCellDirData=data.thetaBasedTimeVsPhaseInfo.(dirStr{di});
            currNumFields=length(fieldnames(currCellDirData));

            currDirTimeBounds=data.cycleTimeFieldBoundPerDirPerField.(dirStr{di});

            for fii=1:currNumFields
                currFieldTimeBound=currDirTimeBounds(fii);



                currFieldNumCyclesPerLap=currCellDirData.(sprintf('field%d',fii)).numThetaCyclesInFieldPerLap;
                goodIdxes=currFieldNumCyclesPerLap>0;

                absoluteCycleTimePerLap=currFieldNumCyclesPerLap(goodIdxes);
                 normalizedCycleTimePerLap=currFieldNumCyclesPerLap(goodIdxes)/currFieldTimeBound;

                 allAbsoluteFieldDurationsPerLap=[allAbsoluteFieldDurationsPerLap;absoluteCycleTimePerLap(:)];
                allNormFieldDurationsPerLap=[allNormFieldDurationsPerLap; normalizedCycleTimePerLap(:)];
            end

        end

    end
    save('durationStatsOct14.mat')
else
    load('durationStatsOct14.mat')
end

onlyFieldBounds=1;
%onlyFieldBounds=0;
minBound=3;

numFields=sum(~isnan(allCycleTimeBounds));

setTightSubplots_Spacious

    useAbsDurations=1;
    %useAbsDurations=0;

    base=1.2;
    base=sqrt(2);


    base=2;
    base=1.333;
        base=1.4;
         %base=1.33333;
     %base=1.1;

    if(useAbsDurations)
        if(onlyFieldBounds)
            histDurations=allCycleTimeBounds(allCycleTimeBounds>minBound);
        else
            histDurations=allAbsoluteFieldDurationsPerLap;
        end
        
         if(onlyFieldBounds)
          maxDur=25;
           maxDur=27;
           %maxDur=50;
          numBins=30;
          %numBins=50;
         else
             maxDur=30;
              numBins=200;
         end
          %maxDur=30;
          %maxDur=40;
        durationEdges=linspace(0,maxDur,numBins);
        
        if(onlyFieldBounds)
            logMin=0;
            logMax=1.2;
        else
            logMin=-1.2;
            logMax=1.2;
        end
        logDurationEdges=linspace(logMin,logMax,numBins);
    else

        histDurations=allNormFieldDurationsPerLap;
          maxDur=5;
        durationEdges=linspace(0,maxDur,50);
        logDurationEdges=linspace(-1,1,100);

    end

    figure;


    durationBins=edgesToBins(durationEdges);

    logDurationBins=edgesToBins(logDurationEdges);

    %histDurations(histDurations>maxDur)=NaN;
    noRectify=1;
    noLimit=1;
    logDurations=getLogTime(histDurations/maxDur,base,noRectify,noLimit);

    Nc=histcounts(histDurations,durationEdges);
    
    
    goodIdxes=(Nc~=0);


    %Nc(Nc==0)=NaN;
    subplot(2,2,1)
    pD=smooth(Nc)/sum(smooth(Nc));
    %plot(durationBins(goodIdxes),pD(goodIdxes),'k-','LineWidth',3)
    plot(durationBins,pD,'k-','LineWidth',4)
xlabel('Field duration (no. of theta cycles)')
ylabel('Probability')
    hold on
box off
    %histfit(histDurations)
    %histogram(allNormFieldDurationsPerLap,linspace(0,7,100))
    axis tight
    %set(gca,'xscale','log')

%{
    subplot(2,1,2)
    %cdfplot(histDurations)
    cdfNc=cumsum(pD);


    plot(durationBins,cdfNc,'k-','LineWidth',3)
    set(gca,'yscale','log')
    set(gca,'xscale','log')
    xlim([0 maxDur])
    %}
    title(sprintf('Field duration distribution (n=%d fields)',numFields))
        subplot(2,2,2)
        
    normplot(histDurations)
    xlabel('Field duration (no. of theta cycles)')
    %box off
    subplot(2,2,3)
    Ncl=histcounts(logDurations,logDurationEdges);
    pcl=smooth(Ncl,3)/sum(smooth(Ncl,3));
    goodIdxes=(Ncl~=0);
    %plot(logDurationBins(goodIdxes),pcl(goodIdxes),'k-','LineWidth',3)
        plot(logDurationBins,pcl,'k-','LineWidth',4)
box off
    %histogram(log(allNormFieldDurationsPerLap),linspace(-1.5,6,100))
    hold on

    %pd = makedist('Normal','mu',0.245,'sigma',0.195)
     pd = makedist('Normal','mu',0.58,'sigma',0.17)
    %pd = makedist('Normal','mu',0.35,'sigma',0.15)
    %pd = makedist('Normal','mu',0.7,'sigma',0.08)

    normFit=pdf(pd,logDurationBins);
      xlabel('log_{1.33}(field duration/max duration)')
    ylabel('Probability')
    
    %yyaxis right
    hold on
    %pL=plot(logDurationBins,normFit/sum(normFit),'r-','LineWidth',2)
    %}
    axis tight
    %legend(pL,'gaussian')
    %legend boxoff
          title(sprintf('log_{%.1f}(duration) distribution (n=%d fields)',base,numFields))

      subplot(2,2,4)
    normplot(logDurations)
        xlabel(sprintf('log_{%.1f}(field duration)',base))
        
        setFigFontTo(16)
        %box off
        saveas(gcf,'fieldDurationLogNormalDistribution.png')
    %{
    yyaxis right
    histfit(logDurations)
    %}







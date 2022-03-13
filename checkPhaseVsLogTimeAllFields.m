%another property of log compression is linear relationship between time and width 
close all; clc

clearvars -except data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load time and phase data within good field traversals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
setTightSubplots_Spacious

numTimeBins=50;
numTimeBins=20;
numTimeBins=15;

numTimeBins=25;
numTimeBins=30;
%numTimeBins=50;
%numTimeBins=35;

timeMin=0.025;
%timeMin=0;

timeMin=0.15;
timeMin=0.05;
timeMin=0.1;
timeMin=0.01;
timeMin=0;
if(timeMin==0)
    timeMinLog=0;
else
    timeMinLog=getLogTime(timeMin);
end

timeFieldEdges=linspace(timeMin,1,numTimeBins+1);
logTimeFieldEdges=linspace(timeMinLog,1,numTimeBins+1);

numPhaseBins=50;
numPhaseBins=120;
numPhaseBins=60;
%numPhaseBins=30;
%numPhaseBins=90;
phaseEdges=linspace(0,360,numPhaseBins+1);


if(~exist('Feb9_2022_phaseVsTimeAndLogTimeStats.mat','file'))
    
    %data=load('allUnitInFieldGoodSpeedSpikePhasesAndTimes.mat');
    tic
    disp('loading spike data....')
    if(~exist('data','var'))
        data=load('allUnitInFieldGoodSpeedSpikePhasesAndTimes.mat');
    end
	toc

    saveFieldPhaseVsLogTimeDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3Images/fieldPhaseVsLogTimeDir';
    touchDir(saveFieldPhaseVsLogTimeDir)

    saveFigureDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/figurePanels';
    touchDir(saveFigureDir)

    %totalUsedFieldCount=data.totalUsedFieldCount;
    totalFieldCount=data.totalFieldCount;

    allFieldSpikePhases=data.allUnitSpikePhases;
    
    %spikeTimeFracInFieldPerField=data.allUnitSpikeTimesInField;
    %{
    allRenormedSpikeTimeFracs=[];
    for fi=1:totalFieldCount
        currFieldRenormedSpikes=data.autoRenormedSpikeTimeFracInFieldPerField{fi};
        
        allRenormedSpikeTimeFracs=[allRenormedSpikeTimeFracs; currFieldRenormedSpikes(:)];
    end
    spikeTimeFracInFieldPerField=allRenormedSpikeTimeFracs;
    %}
    
    allSpikeTimeFracs=[];
    for fi=1:totalFieldCount
        currFieldSpikes=data.spikeTimeFracInFieldPerField{fi};
        
        allSpikeTimeFracs=[allSpikeTimeFracs; currFieldSpikes(:)];
    end
    spikeTimeFracInFieldPerField=allSpikeTimeFracs;

    showPlots=1;
    
        %allFieldToFieldEndTimes=1-spikeTimeFracInFieldPerField;
        
        allFieldToFieldEndTimes=spikeTimeFracInFieldPerField;

        
        %domain restriction also in logTime, to make fair comparison!
        allFieldToFieldEndTimes(allFieldToFieldEndTimes<=0)=NaN;
         allFieldToFieldEndTimes(allFieldToFieldEndTimes>1)=NaN;

         base=1.475;
         base=sqrt(2);
          base=1.33333333;
          
         base=1.4;
         
         base=1.45;
        %allFieldToFieldEndLogTimes=getLogTime(1-spikeTimeFracInFieldPerField-0.15,base);
        timeOffset=0.125;
         timeOffset=0.05;
         timeOffset=0;
        allFieldToFieldEndLogTimes=1-getLogTime(1-spikeTimeFracInFieldPerField-timeOffset,base);
        %allFieldToFieldEndLogTimes=getLogTime(spikeTimeFracInFieldPerField,base);



        if(showPlots)
            figure;
            title(sprintf('Phase vs log(time to end of field), all fields'))
            xlabel('log(time to end of field)')
            ylabel('Theta phase (deg)')
        end

        notNaNIdxes=~isnan(allFieldToFieldEndLogTimes) & ~isnan(allFieldSpikePhases) & timeMinLog<=allFieldToFieldEndLogTimes;

            %[m,b,R]=getLinearFit(currFieldToFieldEndLogTimes(notNaNIdxes),currFieldSpikePhases(notNaNIdxes));
            
            targetSlopeCyclesPerField=1;
            targetOffset=0;
            
            %[ RLogTime,p,s,b ] = kempter_lincirc(currFieldToFieldEndLogTimes(notNaNIdxes),2*pi/360*currFieldSpikePhases(notNaNIdxes),targetSlopeCyclesPerField,targetOffset);
            [ RLogTime,p,s,b ] = kempter_lincirc(allFieldToFieldEndLogTimes(notNaNIdxes),2*pi/360*allFieldSpikePhases(notNaNIdxes));


            mLogTime=s*360; %cycles/field --> degrees/field
            bLogTime=mod(b*360/(2*pi),360); %offset, convert radians to degrees

             notNaNIdxes=~isnan(allFieldToFieldEndTimes) & ~isnan(allFieldSpikePhases) & timeMin<=allFieldToFieldEndTimes;
           %[ RBehavTime,p,s,b ] = kempter_lincirc(currFieldToFieldEndTimes(notNaNIdxes),2*pi/360*currFieldSpikePhases(notNaNIdxes),targetSlopeCyclesPerField,targetOffset);
          [ RBehavTime,p,s,b ] = kempter_lincirc(allFieldToFieldEndTimes(notNaNIdxes),2*pi/360*allFieldSpikePhases(notNaNIdxes));


            mBehavTime=s*360; %cycles/field --> degrees/field
            bBehavTime=mod(b*360/(2*pi),360); %offset, convert radians to degrees

         

        if(showPlots)
            fH=figure; 
            subplot(1,2,2)
            %plot(allFieldToFieldEndLogTimes,allFieldSpikePhases,'k.')

            %getJointDistrGivenX(allFieldToFieldEndLogTimes,allFieldSpikePhases,logTimeFieldEdges,phaseEdges,fH);
            getJointDistr(allFieldToFieldEndLogTimes,allFieldSpikePhases,logTimeFieldEdges,phaseEdges,fH);

            xlabel('log(time to end of field) (normalized)')
             %xlabel('log(1-time in field)')
            ylabel('Theta phase (deg)')

            %xlim([0 360])
            %xlim([-20 380])

            %title(sprintf('Theta phase  vs log(time to end), field %d (n=%d spikes)',fi,length(currFieldSpikePhases)))


            %daspect([1 360 1])
            %caxis([0 prctile(pxySmooth(:),97.5)])
            %maxFig
            %subplot(1,2,2)


            hold on
            x0=0;
            y0=bLogTime;
            xf=1;
            yf=mLogTime+bLogTime;

            yMid=(y0+yf)/2;



            %dataMid=circMedianDeg(allFieldSpikePhases);


            bestShift=NaN;
            shifts=-5:5;
            yMids=NaN(size(shifts));
            y0s=NaN(size(shifts));
            yFs=NaN(size(shifts));
            shiftErrors=NaN(size(shifts));
            for si=1:length(shifts)

                y0Temp=y0+360*shifts(si);
               yfTemp=yf+360*shifts(si);

               %for i=1:length(currFieldSpikePhases)
               predictedVals=interp1([0 1],[y0Temp yfTemp],allFieldToFieldEndLogTimes);
               actualVals=allFieldSpikePhases;

               errors=abs(actualVals-predictedVals);


               meanError=nanmean(errors);

               shiftErrors(si)=meanError;
               %end


               %{
                yMids(si)=(y0Temp+yfTemp)/2;
                yFs(si)=yfTemp;
                y0s(si)=y0Temp;
                %}
            end

             errors=rad2ang(angdiff(ang2rad(actualVals),ang2rad(predictedVals)));
           
             %meanErrorLogTime=nanmean(abs(errors));
             meanErrorLogTime=nanmean((errors));
             %meanErrorLogTime=circMeanDeg(errors);

            %figure; histogram(errors,30)

            %[~,bestShift]=min(abs(dataMid-yMids) + abs(360-yFs)+ abs(0-y0s));
           [~,bestShift]=min(abs(shiftErrors));


            y0=y0+360*shifts(bestShift);
               yf=yf+360*shifts(bestShift);


            plot([x0 xf], [y0 yf],'k--','LineWidth',5)
            plot([x0 xf], [y0 yf]+360,'k--','LineWidth',5)
                 plot([x0 xf], [y0 yf]-360,'k--','LineWidth',5)
                 
            title({sprintf('Theta phase vs log(time to end), all fields (n=%d spikes)',length(allFieldSpikePhases)),sprintf('circR=%.2f, slope=%.2f deg/field, mean circ error = %.2f deg',RLogTime,mLogTime,meanErrorLogTime)})

            
             %xlabel('log(1- (time in field))')
             xlabel('log(time in field)')
            ylabel('Theta phase (deg)')

            ylim([0 360])
            xlim([timeMinLog 1])
            axis square
            %daspect([0 360 1])
             box off
             
             subplot(1,2,1)
            %plot(allFieldToFieldEndTimes,allFieldSpikePhases,'k.')
            %getJointDistrGivenX(allFieldToFieldEndTimes,allFieldSpikePhases,timeFieldEdges,phaseEdges,fH);
            getJointDistr(allFieldToFieldEndTimes,allFieldSpikePhases,timeFieldEdges,phaseEdges,fH);

             xlabel('Time in field (normalized)')
            ylabel('Theta phase (deg)')


             x0=0;
            y0=bBehavTime;
            xf=1;
            yf=mBehavTime+bBehavTime;

             bestShift=NaN;
            shifts=-5:5;
            yMids=NaN(size(shifts));
            y0s=NaN(size(shifts));
            yFs=NaN(size(shifts));
            shiftErrors=NaN(size(shifts));
            
            for si=1:length(shifts)
                y0Temp=y0+360*shifts(si);
               yfTemp=yf+360*shifts(si);

               %for i=1:length(currFieldSpikePhases)
               predictedVals=interp1([0 1],[y0Temp yfTemp],allFieldToFieldEndLogTimes);
               actualVals=allFieldSpikePhases;

               errors=abs(actualVals-predictedVals);
               %errors=rad2ang(angdiff(ang2rad(actualVals),ang2rad(predictedVals)));

               meanError=nanmean(errors);

               shiftErrors(si)=meanError;
               %end

            end

            %[~,bestShift]=min(abs(dataMid-yMids) + abs(360-yFs)+ abs(0-y0s));
           [~,bestShift]=min(abs(shiftErrors));

            errors=rad2ang(angdiff(ang2rad(actualVals),ang2rad(predictedVals)));



                    %meanErrorBehavTime=nanmean(abs(errors));
                    meanErrorBehavTime=nanmean((errors));
               %meanErrorBehavTime=circMeanDeg(errors);



            y0=y0+360*shifts(bestShift);
               yf=yf+360*shifts(bestShift);

            hold on

             plot([x0 xf], [y0 yf],'k--','LineWidth',5)
                 plot([x0 xf], [y0 yf]+360,'k--','LineWidth',5)
                 plot([x0 xf], [y0 yf]-360,'k--','LineWidth',5)

             ylim([0 360])
            xlim([timeMin 1])
             box off
             axis square

            title({sprintf('Theta phase vs Time in field, all fields (n=%d spikes)',length(allFieldSpikePhases)),sprintf('circR=%.2f, slope=%.2f deg/field, mean circ error = %.2f deg',RBehavTime,mBehavTime,meanErrorBehavTime)})


            %set(gca,'xscale','log')
            setFigFontTo(18)
            %maxFigHalfWidth
            maxFig

            saveas(gcf,fullfile(saveFieldPhaseVsLogTimeDir,sprintf('logTimeVsThetaPhase_AllFields.png')))
            close all

        end

       


    

    %%

    %save('phaseVsTimeAndLogTimeStats.mat')

else
    load('Feb9_2022_phaseVsTimeAndLogTimeStats.mat')
    %load('phaseVsTimeAndLogTimeStats.mat')
end
    
%%
%numFieldsUsed
setTightSubplots_Spacious

numRhoBins=50;
rhoEdges=linspace(-1,1,numRhoBins+1);
rhoBins=edgesToBins(rhoEdges);

numErrorBins=50;
numErrorBins=100;
minError=-90;
maxError=90;

minError=-80;
maxError=80;
errorEdges=linspace(minError,maxError,numErrorBins+1);

errorBins=edgesToBins(errorEdges);

figure
subplot(2,1,1)
Nlog=histcounts(allFieldLogTimePhaseRs,rhoEdges);
Nlin=histcounts(allFieldTimePhaseRs,rhoEdges);
pLinRho=plot(rhoBins,smooth(Nlin),'k-','LineWidth',4)
hold on
pLogRho=plot(rhoBins,smooth(Nlog),'r-','LineWidth',4)

title(sprintf('circR distributions, %d fields',numFieldsUsed ))
xlim([-1 1])

legend([pLinRho pLogRho],{'phase-time linear model','phase-log(time) linear model'},'Location','northwest')

xlabel('Circular correlation coefficient')
ylabel('Field count')

box off
subplot(2,1,2)
NlogErr=histcounts(allFieldLogTimeMeanCircErrors,errorEdges);
NlinErr=histcounts(allFieldTimeMeanCircErrors,errorEdges);
pLinErr=plot(errorBins,smooth(NlinErr),'k-','LineWidth',4)
hold on
pLogErr=plot(errorBins,smooth(NlogErr),'r-','LineWidth',4)
xlabel('Mean error (deg)')
ylabel('Field count')
legend([pLinErr pLogErr],{'phase-time linear model','phase-log(time) linear model'},'Location','northwest')
box off

title(sprintf('error distributions, %d fields',numFieldsUsed ))

setFigFontTo(18)
maxFigHalfWidth

saveas(gcf,'linearModelLogTimeFitStatsDistr.png')
%xlim([-1 1])

%{
histogram(allFieldLogTimePhaseSlopes,30)
title(sprintf('circSlope distribution (deg/field), %d fields',numFieldsUsed ))

xlim([-360 720])
%}






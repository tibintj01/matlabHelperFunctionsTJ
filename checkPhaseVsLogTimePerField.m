%log fit of phase vs time, appropriate time offset and log base
clear all; close all; clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load time and phase data within good field traversals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
setTightSubplots_Spacious

useSummaryPhasePerCycle=1;

threshBehavR=0.8;
threshBehavR=0.7;
%threshBehavR=0.5;
threshBehavR=0;

%threshBehavSlope=270; %deg per field
%maxBehavSlope=450;
%threshBehavSlope=180; %deg per field
%threshBehavSlope=240; %deg per field
%maxBehavSlope=480;
threshBehavSlope=270; %deg per field
maxBehavSlope=450;

slopeRange=180;
%slopeRange=360;
slopeRange=90;
threshBehavSlope=360-slopeRange; %deg per field
maxBehavSlope=360+slopeRange;


numFieldsUsed=0;
%{
threshBehavSlope=260; %deg per field
maxBehavSlope=460;
%}

%{
threshBehavSlope=280; %deg per field
maxBehavSlope=440;

threshBehavSlope=288; %deg per field
maxBehavSlope=432;
%}

%threshBehavSlope=300; %deg per field
%maxBehavSlope=420;

%threshBehavSlope=210; %deg per field
%maxBehavSlope=510;
%threshBehavSlope=180; %deg per field
%maxBehavSlope=540;

    showPlots=0;
    %showPlots=1;
    
    fD=figure(191);
    
  sourceDataFilePath='allUnitInFieldGoodSpeedSpikePhasesAndTimes.mat';

if(~exist('dphaseVsTimeAndLogTimeStats.mat','file'))
    
    data=load(sourceDataFilePath)

    saveFieldPhaseVsLogTimeDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3Images/fieldPhaseVsLogTimeDir';
    touchDir(saveFieldPhaseVsLogTimeDir)

    saveFigureDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/figurePanels';
    touchDir(saveFigureDir)

    
    %trueTimeBoundFrac=0.8;
    trueTimeBoundFrac=1;
    %totalUsedFieldCount=data.totalUsedFieldCount;
    totalFieldCount=data.totalFieldCount;

    spikePhasesPerField=data.spikePhasesPerField;
    spikeTimeFracInFieldPerField=data.spikeTimeFracInFieldPerField;
        spikeTimeFracInFieldPerField=data.spikeTimeFracInFieldPerFieldRenormed;

        
    %spikeTimeFracInFieldPerFieldRenormed=spikeTimeFracInFieldPerField;
    
        %spikeTimeFracInFieldPerField=data.autoRenormedSpikeTimeFracInFieldPerField;


    minNumSpikesForField=600;
    minNumSpikesForField=700;
    minNumSpikesForField=300;

    minNumSpikesForField=30;
       % minNumSpikesForField=50;
    %minNumSpikesForField=3000;

    numFieldsUsed=0;

    allFieldLogTimePhaseSlopes=[];
    allFieldLogTimePhaseRs=[];

    allInvPhaseGainR=[];
    allInvPhaseGainSlopes=[];

    allFieldTimePhaseSlopes=[];
    allFieldTimePhaseRs=[];

    allFieldLogTimeMeanCircErrors=[];
    allFieldTimeMeanCircErrors=[];
    
    allFieldDerivWithTimeVals=[];
    allFieldDerivWithTimePts=[];
    
    allFieldsLogErrVsTime=[];
    allFieldsLinErrVsTime=[];
    
    allFieldsLogCentralConcOfErr=[];
    allFieldsLinCentralConcOfErr=[];

    
    %xDispMin=-1;
     xDispMin=0;
     
     fieldIsUsedForLogTimeStats=zeros(totalFieldCount,1);

     fiStart=1;
    for fi=fiStart:totalFieldCount

        currFieldSpikePhases=spikePhasesPerField{fi};
        currFieldInFieldTimes=spikeTimeFracInFieldPerField{fi};
        if(size(currFieldInFieldTimes,2)>1)
            currFieldInFieldTimes=currFieldInFieldTimes(:,1);
        end
        %spikeTimeFracInFieldPerFieldRenormed{fi}

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %RESCALE MANUALLY CHOSEN TIME BOUNDS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        currFieldInFieldTimes=currFieldInFieldTimes/trueTimeBoundFrac;
        
        currFieldToFieldEndTimes=1-currFieldInFieldTimes;
         %currFieldToFieldEndTimes=currFieldInFieldTimes;
        %currFieldToFieldEndTimes=1-currFieldInFieldTimes;

        
        %domain restriction also in logTime, to make fair comparison!
        %currFieldToFieldEndTimes(currFieldToFieldEndTimes<=xDispMin)=NaN;
        currFieldToFieldEndTimes(currFieldToFieldEndTimes<0)=NaN;

         currFieldToFieldEndTimes(currFieldToFieldEndTimes>1)=NaN;

         
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %choosing offset (location of effective zero) and log base?
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
         base=2;
         
        base=sqrt(2);
        base=sqrt(2);
           %base=2;
        base=1.7;
        base=1.335;
         base=1.31;
         base=1.3333333;
         base=1.4;
         
        % base=1.5457;
             %base=sqrt(2);
          %base=1.34;
         %timeOffset=0.05;
          %timeOffset=-0.05;
           timeOffset=0;
           %timeOffset=-0.05;
            %timeOffset=0.1;
        currFieldToFieldEndLogTimes=getLogTime(1-currFieldInFieldTimes-timeOffset,base);

        if(length(currFieldSpikePhases)<minNumSpikesForField)
            continue
        else
            %numFieldsUsed=numFieldsUsed+1;
            %continue
        end


       

        notNaNIdxes=~isnan(currFieldToFieldEndLogTimes) & ~isnan(currFieldSpikePhases) & currFieldToFieldEndLogTimes>=xDispMin;

        
        if(useSummaryPhasePerCycle)
            currFieldLogTimeCenters=unique(currFieldToFieldEndLogTimes(notNaNIdxes));
            
            currFieldBehavTimeCenters=unique(currFieldToFieldEndTimes(notNaNIdxes));
            
            [summaryPhasePerCycleLogTime] = getSummaryPhasePerCycle(currFieldLogTimeCenters,currFieldToFieldEndLogTimes,currFieldSpikePhases);
            
            [summaryPhasePerCycleBehavTime] = getSummaryPhasePerCycle(currFieldBehavTimeCenters,currFieldToFieldEndTimes,currFieldSpikePhases);

            disp('')
            
        end


            %[m,b,R]=getLinearFit(currFieldToFieldEndLogTimes(notNaNIdxes),currFieldSpikePhases(notNaNIdxes));
            
            targetSlopeCyclesPerField=1;
            targetOffset=0;
            
            %[ RLogTime,p,s,b ] = kempter_lincirc(currFieldToFieldEndLogTimes(notNaNIdxes),2*pi/360*currFieldSpikePhases(notNaNIdxes),targetSlopeCyclesPerField,targetOffset);
             if(useSummaryPhasePerCycle)
                 
                  notNaNIdxes=~isnan(summaryPhasePerCycleLogTime) & ~isnan(currFieldLogTimeCenters);

                %[ RLogTime,p,s,b ] = kempter_lincirc(currFieldLogTimeCenters(notNaNIdxes),2*pi/360*summaryPhasePerCycleLogTime(notNaNIdxes));
                
                 
                    %[bestShiftDeg,shiftedPhases]=getBestShiftDeg1D(summaryPhasePerCycleLogTime(notNaNIdxes),1-currFieldLogTimeCenters(notNaNIdxes),0.3);
                                        [bestShiftDeg,shiftedPhases]=getBestShiftDeg1D(summaryPhasePerCycleLogTime(notNaNIdxes),1-currFieldLogTimeCenters(notNaNIdxes),0.3);

                 currFieldSpikePhases=mod(currFieldSpikePhases-bestShiftDeg,360);
                 summaryPhasePerCycleLogTime=mod(summaryPhasePerCycleLogTime-bestShiftDeg,360);
                 
                    [mLogTime,bLogTime,RLogTime]=getLinearFit(currFieldLogTimeCenters(notNaNIdxes),summaryPhasePerCycleLogTime(notNaNIdxes));
                    
             else
                [ RLogTime,p,s,b ] = kempter_lincirc(currFieldToFieldEndLogTimes(notNaNIdxes),2*pi/360*currFieldSpikePhases(notNaNIdxes));
                  mLogTime=s*360; %cycles/field --> degrees/field
                bLogTime=mod(b*360/(2*pi),360); %offset, convert radians to degrees
             end

            
          
            

             notNaNIdxes=~isnan(currFieldToFieldEndTimes) & ~isnan(currFieldSpikePhases) & currFieldToFieldEndTimes>=xDispMin;
           %[ RBehavTime,p,s,b ] = kempter_lincirc(currFieldToFieldEndTimes(notNaNIdxes),2*pi/360*currFieldSpikePhases(notNaNIdxes),targetSlopeCyclesPerField,targetOffset);
          if(useSummaryPhasePerCycle)
              notNaNIdxes=~isnan(summaryPhasePerCycleBehavTime) & ~isnan(currFieldBehavTimeCenters);
              
            %[ RBehavTime,p,s,b ] = kempter_lincirc(currFieldBehavTimeCenters(notNaNIdxes),2*pi/360*summaryPhasePerCycleBehavTime(notNaNIdxes));
             %[bestShiftDeg,shiftedPhases]=getBestShiftDeg1D(summaryPhasePerCycleBehavTime(notNaNIdxes),1-currFieldBehavTimeCenters(notNaNIdxes),0.3);
                 %currFieldSpikePhases=mod(currFieldSpikePhases-bestShiftDeg,360);
                 
                 summaryPhasePerCycleBehavTime=mod(summaryPhasePerCycleBehavTime-bestShiftDeg,360);
                 
            [mBehavTime,bBehavTime,RBehavTime]=getLinearFit(currFieldBehavTimeCenters(notNaNIdxes),summaryPhasePerCycleBehavTime(notNaNIdxes));
            

          else
           [ RBehavTime,p,s,b ] = kempter_lincirc(currFieldToFieldEndTimes(notNaNIdxes),2*pi/360*currFieldSpikePhases(notNaNIdxes));
            mBehavTime=s*360; %cycles/field --> degrees/field
            bBehavTime=mod(b*360/(2*pi),360); %offset, convert radians to degrees
          end

     
          if(RBehavTime<threshBehavR || mBehavTime<threshBehavSlope || mBehavTime>maxBehavSlope)
              continue 
          end

            if(isnan(RLogTime))
                continue
            end


           

        if(showPlots)
            
            figure;
            title(sprintf('Phase vs log(time to end of field), field %d',fi))
            xlabel('log(time to end of field)')
            ylabel('Theta phase (deg)')
  
                subplot(1,2,1)
                plot(1-currFieldToFieldEndLogTimes,currFieldSpikePhases,'k.','MarkerSize',15)

                if(useSummaryPhasePerCycle)
                    hold on

                    plot(1-currFieldLogTimeCenters,summaryPhasePerCycleLogTime,'ro','MarkerSize',15)
                     plot(1-currFieldLogTimeCenters,summaryPhasePerCycleLogTime,'r.','MarkerSize',20)
                end

                 xlabel('log(time to end of field) (normalized)')
                ylabel('Theta phase (deg)')

                %xlim([0 360])
                %xlim([-20 380])

                %title(sprintf('Theta phase  vs log(time to end), field %d (n=%d spikes)',fi,length(currFieldSpikePhases)))


                %daspect([1 360 1])
                %caxis([0 prctile(pxySmooth(:),97.5)])
                %maxFig
                %subplot(1,2,2)


                hold on
                
                
           end
        
            x0=xDispMin;
            y0=mLogTime*xDispMin+bLogTime;
            xf=1;
            yf=mLogTime+bLogTime;

            yMid=(y0+yf)/2;


            dataMid=circMedianDeg(currFieldSpikePhases);


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
               predictedVals=interp1([0 1],[y0Temp yfTemp],currFieldToFieldEndLogTimes);
               actualVals=currFieldSpikePhases;

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
             %if(isempty(errors))
             %    continue
             %end
             %meanErrorLogTime=nanmean(abs(errors));
             meanErrorLogTime=nanmean((errors));
             
             currFieldLogErrVsTime=[currFieldToFieldEndLogTimes(:) errors(:)];
             
             currFieldCentralConcOfErrLog=getCentralConcentrationOfErr(invertLogTime(currFieldToFieldEndLogTimes(:),base),errors(:));
             
             allFieldsLogCentralConcOfErr=[allFieldsLogCentralConcOfErr; currFieldCentralConcOfErrLog];
             
             allFieldsLogErrVsTime=[allFieldsLogErrVsTime; currFieldLogErrVsTime];
             %meanErrorLogTime=circMeanDeg(errors);

            %figure; histogram(errors,30)

            %[~,bestShift]=min(abs(dataMid-yMids) + abs(360-yFs)+ abs(0-y0s));
           [~,bestShift]=min(abs(shiftErrors));


            y0=y0+360*shifts(bestShift);
               yf=yf+360*shifts(bestShift);

              if(showPlots)

                    plot(1-[x0 xf], [y0 yf],'k--','LineWidth',5)
                    plot(1-[x0 xf], [y0 yf]+360,'k--','LineWidth',5)
                         plot(1-[x0 xf], [y0 yf]-360,'k--','LineWidth',5)

                    title({sprintf('Theta phase vs log(time to end), field %d (n=%d spikes)',fi,length(currFieldSpikePhases)),sprintf('circR=%.2f, slope=%.2f deg/field, mean circ error = %.2f deg',RLogTime,mLogTime,meanErrorLogTime)})

                     xlabel('log(time to end of field) (normalized)')
                    ylabel('Theta phase (deg)')

                    ylim([0 360])
                    xlim([xDispMin 1])
                     box off
                     subplot(1,2,2)
                    plot(1-currFieldToFieldEndTimes,currFieldSpikePhases,'k.','MarkerSize',15)
                    if(useSummaryPhasePerCycle)
                        hold on
                        plot(1-currFieldBehavTimeCenters,summaryPhasePerCycleBehavTime,'ro','MarkerSize',15)

                        plot(1-currFieldBehavTimeCenters,summaryPhasePerCycleBehavTime,'r.','MarkerSize',20)

                    end
                     xlabel('Time to end of field (normalized)')
                    ylabel('Theta phase (deg)')

              end
               
              dt=median(diff(currFieldBehavTimeCenters));
              
              if(dt>0.25) %minimum theta time diff
                  continue
              end
              derivativeEstWind= floor(0.2/dt);%corresponding to approximately 0.2 of field;
               %derivativeEstWind= ceil(0.25/dt);%quad approximation corresponding to 0.4 of field;

              if(derivativeEstWind<=2)
                  derivativeEstWind=3; %minumum window 
              end
             
              derivModelOrder=2; %local quadratic fit
               currFieldPhaseGainsPerTime=movingslope(summaryPhasePerCycleBehavTime,derivativeEstWind,derivModelOrder,dt);
               
               %remove edge effects
               
               currFieldPhaseGainsPerTime(1)=NaN;
               currFieldPhaseGainsPerTime(end)=NaN; 
               
               
               currFieldPhaseGainsPerTimePts=currFieldBehavTimeCenters;
               
               %{
               figure; plot(currFieldBehavTimeCenters,summaryPhasePerCycleBehavTime,'r.','MarkerSize',20)
               yyaxis right;plot(currFieldBehavTimeCenters,currFieldPhaseGainsPerTime,'k.','MarkerSize',30)
               ylim([0 Inf])
               close all
               %}
               
            x0=xDispMin;
            y0=mBehavTime*xDispMin+bBehavTime;
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
               %predictedVals=interp1([0 1],[y0Temp yfTemp],currFieldToFieldEndLogTimes);
               predictedVals=interp1([0 1],[y0Temp yfTemp],currFieldToFieldEndTimes);
               
               actualVals=currFieldSpikePhases;

               errors=abs(actualVals-predictedVals);
               %errors=rad2ang(angdiff(ang2rad(actualVals),ang2rad(predictedVals)));

               meanError=nanmean(errors);


               shiftErrors(si)=meanError;
               %end


               %{
                yMids(si)=(y0Temp+yfTemp)/2;
                yFs(si)=yfTemp;
                y0s(si)=y0Temp;
                %}
            end

            %[~,bestShift]=min(abs(dataMid-yMids) + abs(360-yFs)+ abs(0-y0s));
           [~,bestShift]=min(abs(shiftErrors));

            errors=rad2ang(angdiff(ang2rad(actualVals),ang2rad(predictedVals)));


            currFieldCentralConcOfErrLin=getCentralConcentrationOfErr(currFieldToFieldEndTimes(:),errors(:));
            allFieldsLinCentralConcOfErr=[allFieldsLinCentralConcOfErr; currFieldCentralConcOfErrLin];
            
             %if(isempty(errors))
             %    continue
             %end

                    %meanErrorBehavTime=nanmean(abs(errors));
                    meanErrorBehavTime=nanmean((errors));
               %meanErrorBehavTime=circMeanDeg(errors);

               
               currFieldLinErrVsTime=[currFieldToFieldEndTimes(:) errors(:)];
             
             allFieldsLinErrVsTime=[allFieldsLinErrVsTime; currFieldLinErrVsTime];

             

            y0=y0+360*shifts(bestShift);
               yf=yf+360*shifts(bestShift);
               
                     %only from initial max value of gain
                    [~,maxPhaseGainID]=max(currFieldPhaseGainsPerTime);
                    
                    %{
                    if(maxPhaseGainID>2)
                        currFieldPhaseGainsPerTime(1:(maxPhaseGainID)-2)=NaN;
                    end
                    %}
                    
                    
                      if(maxPhaseGainID>1)
                        currFieldPhaseGainsPerTime(1:(maxPhaseGainID)-1)=NaN;
                    end
                    
                    
                    
                    currFieldPhaseGainsPerTime(currFieldPhaseGainsPerTime<0)=NaN;
                    

               if(showPlots)
                    hold on

                     plot(1-[x0 xf], [y0 yf],'k--','LineWidth',5)
                         plot(1-[x0 xf], [y0 yf]+360,'k--','LineWidth',5)
                         plot(1-[x0 xf], [y0 yf]-360,'k--','LineWidth',5)

                     ylim([0 360])
                    xlim([xDispMin 1])
                     box off

                    title({sprintf('Theta phase vs Time to end, field %d (n=%d spikes)',fi,length(currFieldSpikePhases)),sprintf('circR=%.2f, slope=%.2f deg/field, mean circ error = %.2f deg',RBehavTime,mBehavTime,meanErrorBehavTime)})


                     setFigFontTo(18)
                    %maxFigHalfWidth
                    maxFig
                    
                    %{
                    [newFieldEnd,~] = ginput
                %assert(length(newFieldEnd)==1)
                
                if(isempty(newFieldEnd))
                    newNorm=1;
                else
                    newNorm=1-newFieldEnd;
                end
                spikeTimeFracInFieldPerFieldRenormed{fi}=spikeTimeFracInFieldPerFieldRenormed{fi}/newNorm;
                    %}

                    saveas(gcf,fullfile(saveFieldPhaseVsLogTimeDir,sprintf('logTimeVsThetaPhase_Field%d.png',fi)))
                   
                    
              
                    
                    
                    %{
                     figure
                     subplot(2,1,1)
                    %plot(allFieldDerivWithTimePts,allFieldDerivWithTimeVals,'k.','MarkerSize',15)
                    plot(currFieldPhaseGainsPerTimePts,currFieldPhaseGainsPerTime,'k.','MarkerSize',15)
                    ylim([0 1800])
                    hold on
                       xlim([0 1])
                    subplot(2,1,2)
                    %plot(allFieldDerivWithTimePts,1./allFieldDerivWithTimeVals,'k.','MarkerSize',15)
                    plot(currFieldPhaseGainsPerTimePts,1./currFieldPhaseGainsPerTime,'k.','MarkerSize',15)
                    ylim([0 0.01])
                    
                    xlim([0 1])
                    %}
                    close all
                    
               end
          

  
               derivNoNaNIdxes=~isnan(currFieldPhaseGainsPerTimePts) & ~isnan(1./currFieldPhaseGainsPerTime) & ~isinf(1./currFieldPhaseGainsPerTime);
                y=1./((currFieldPhaseGainsPerTime(derivNoNaNIdxes)/360));
                x=1-currFieldPhaseGainsPerTimePts(derivNoNaNIdxes);

               %allFieldDerivWithTimePts
               %1./allFieldDerivWithTimeVals
               [invPhaseGainSlope,b,invPhaseGainR]=getLinearFit(x,y);
               
               
               
          allInvPhaseGainR=[allInvPhaseGainR;invPhaseGainR];
          allInvPhaseGainSlopes=[allInvPhaseGainSlopes;invPhaseGainSlope];

        allFieldLogTimePhaseSlopes=[allFieldLogTimePhaseSlopes; mLogTime];
        allFieldLogTimePhaseRs=[allFieldLogTimePhaseRs; RLogTime];

        allFieldTimePhaseSlopes=[allFieldTimePhaseSlopes; mBehavTime];
        allFieldTimePhaseRs=[allFieldTimePhaseRs; RBehavTime];


        allFieldLogTimeMeanCircErrors=[allFieldLogTimeMeanCircErrors;meanErrorLogTime];
        allFieldTimeMeanCircErrors=[allFieldTimeMeanCircErrors; meanErrorBehavTime];
        
        
        allFieldDerivWithTimeVals=[allFieldDerivWithTimeVals ;currFieldPhaseGainsPerTime(:)];
        allFieldDerivWithTimePts=[allFieldDerivWithTimePts; currFieldPhaseGainsPerTimePts(:)];
    
        length(allFieldTimeMeanCircErrors)
        
        fieldIsUsedForLogTimeStats(fi)=1;
        
        numFieldsUsed=numFieldsUsed+1;
       

        hold on
    end

    %%
    %save(sourceDataFilePath,'spikeTimeFracInFieldPerFieldRenormed','-append')

    save('phaseVsTimeAndLogTimeStats.mat')

else
    load('phaseVsTimeAndLogTimeStats.mat')
end
    
%%

figure; 
subplot(1,2,1)
plot(1-allFieldDerivWithTimePts,allFieldDerivWithTimeVals/360,'k.','MarkerSize',8)
ylim([0 1500]/360)

hold on

timeAxis=linspace(0.01,1,100);
timeAxis=linspace(0.15,1,100);
%plot(timeAxis,100./(timeAxis)+100,'k--','LineWidth',3)
%plot(timeAxis,150./(timeAxis)+20,'k--','LineWidth',3)
 %plot(timeAxis,200./(timeAxis)-50,'k--','LineWidth',3)
 %plot(timeAxis,175./(timeAxis)-50,'k--','LineWidth',3)
 %plot(timeAxis,175./(timeAxis-0.05)-50,'r--','LineWidth',5)
 %plot(timeAxis,100./(timeAxis-0.05)+50,'r--','LineWidth',5)
%plot(timeAxis,50./(timeAxis-0.15)+50,'r--','LineWidth',5)
%plot(timeAxis,75./(timeAxis-0.15)+75,'r--','LineWidth',5)
%oneOverXh=plot(1-timeAxis,90./(timeAxis-0.15),'r--','LineWidth',5) 
%oneOverXh=plot(1-timeAxis,(1/4)./(timeAxis-0.15),'r--','LineWidth',5) 
%oneOverXh=plot(1-timeAxis,(1/4)./(timeAxis-0.15),'r--','LineWidth',5)
timeAxis=linspace(0,0.85,100);
timeOffset=0.1;
%oneOverXh=plot(1-timeAxis,(1/2.2964)./(timeAxis),'r--','LineWidth',5)
%oneOverXh=plot(timeAxis,(0.3466)./(1-timeAxis),'r--','LineWidth',5)

%oneOverXh=plot(timeAxis,(log(1.5457))./(1-timeAxis),'r--','LineWidth',5)
%oneOverXh=plot(timeAxis,(log(1.333333))./(1-timeAxis-timeOffset),'r--','LineWidth',5)
oneOverXh=plot(timeAxis,(log(base))./(1-timeAxis-timeOffset),'r--','LineWidth',5)

%exp(1/2.2964)=1.5457

%oneOverXh=plot(1-timeAxis,(1/4)./(timeAxis),'r--','LineWidth',5)
%oneOverXh=plot(1-timeAxis,87.5274./(timeAxis-0.15),'r--','LineWidth',5) 
%oneOverXh=plot(1-timeAxis,87.5274./(timeAxis-0.2),'r--','LineWidth',5) 
%oneOverXh=plot(1-timeAxis,75./(timeAxis-0.2)+25,'r--','LineWidth',5) 
%oneOverXh=plot(1-timeAxis,60./(timeAxis)+75,'r--','LineWidth',5) 

%xlabel('Time to end of field (frac)')
xlabel('Time in field (frac)')
ylabel('Local phase derivative magnitude (cycles/field)')
ylim([0 4])

legend(oneOverXh,'d ~ 1/(1-t)','Location','northwest')

legend boxoff
maxFigHalfWidth
box off

title(sprintf('Local phase gain across field (n=%d fields)',length(allFieldTimeMeanCircErrors)))
xlim([0 1-timeOffset])


%allFieldDerivWithTimePts=allFieldDerivWithTimePts-0.15;
%plot(timeAxis,200./(timeAxis-0.05)-100,'r-','LineWidth',3)
%figure
subplot(1,2,2)

%ylim([0 5])
derivNoNaNIdxes=~isnan(allFieldDerivWithTimePts) & ~isnan(1./allFieldDerivWithTimeVals) & ~isinf(1./allFieldDerivWithTimeVals);
y=1./((allFieldDerivWithTimeVals(derivNoNaNIdxes)/360));
x=1-allFieldDerivWithTimePts(derivNoNaNIdxes);

%plot(allFieldDerivWithTimePts,1./(1-(allFieldDerivWithTimeVals/360)),'k.')
plot(x,y,'k.','MarkerSize',8)
box off


ylim([0 4])
%X=[ones(length(x),1), x(:)];
X=[ x(:)];
xlabel('Time in field (frac)')
ylabel('1/(phase derivative magnitude)')

[robustFitCoeffs,stats] = robustfit(X,y);
hold on

%plot([0 1],robustFitCoeffs,'r--','LineWidth',5)
%plot([0 1],[1/log(1.333333) 0],'r--','LineWidth',5)

plot([0 1]-timeOffset,[1/log(base) 0],'r--','LineWidth',5)



%xlim([timeOffset 1])
%[m,b,R]=getLinearFit(1-X,y,1,0)
xlim([0 1-timeOffset])

%plot(1-[timeOffset 1],[0 1/log(1.333333)],'b--','LineWidth',5)
%plot(1-[0 1],[0 4],'r--','LineWidth',5)
setFigFontTo(14)

saveas(gcf,'phaseSlopeOverFieldOneOverX.png')


%%
 figure
 %subplot(3,1,1)
  rEdges=linspace(-1,1,61);
 rBins=edgesToBins(rEdges);
 Nc=histcounts(allInvPhaseGainR,rEdges); 
 plot(rBins,smooth(Nc,3)/sum(smooth(Nc,3)),'r-','LineWidth',5)
 axis tight
 hold on
 plot([0 0], ylim,'k--','LineWidth',3)
 
xlabel('Inverse phase gain-time R')
ylabel('Probability')
box off

%%%%%%%%%%%%%%%
%STATS
%%%%%%%%%%%%%%%
nanmean(allInvPhaseGainR)
getSEMacrossRows(allInvPhaseGainR)
[p,h,stats]=signrank(allInvPhaseGainR)

figure
ylim([0 35])
 subplot(3,1,2)
 %
  baseEdges=linspace(-10,1,61);
 baseBins=edgesToBins(baseEdges);
 %
  Nc=histcounts(allInvPhaseGainSlopes,baseEdges); 
 plot(baseBins,smooth(Nc)/sum(smooth(Nc)),'r-','LineWidth',5)
 xlabel('Inverse phase gain-time slope (fields/cycle)')
 

 %%%%%%%%%%%%%%%
%STATS
%%%%%%%%%%%%%%%
nanmean(allInvPhaseGainSlopes)
getSEMacrossRows(allInvPhaseGainSlopes)
[p,h,stats]=signrank(allInvPhaseGainSlopes)

    [~,maxID]=max(smooth(Nc));
 bestBase=baseBins(maxID);
 hold on
  plot([bestBase bestBase], ylim,'k--','LineWidth',3)
  
  axis tight

 
 title(sprintf('Maximum probabillity slope: %.2f',bestBase))
 
%ylabel('Field count')
ylabel('Probability')

 subplot(3,1,3)
 
 baseEdges=linspace(0,4,51);
  %baseEdges=linspace(0.5,3,41);
  baseBins=edgesToBins(baseEdges);
 Nc=histcounts(exp(-1./allInvPhaseGainSlopes),baseEdges); 
  plot(baseBins,smooth(Nc,3)/sum(smooth(Nc,3)),'r-','LineWidth',5)
  
       [~,maxID]=max(smooth(Nc))
 %bestBase=exp(-1/baseBins(maxID))
  %bestBase=baseBins(maxID);
 axis tight
  hold on
 plot(exp(-1./[bestBase bestBase]), ylim,'k--','LineWidth',3)



title(sprintf('Maximum probabillity logarithm base: %.3f',exp(-1/bestBase)))
xlabel('Best fit logarithm base')

%ylabel('Field count')
ylabel('Probability')
setFigFontTo(14)
box off
axis tight
numFieldsUsed

%%%%%%%%%%%%%%%
%STATS
%%%%%%%%%%%%%%%
nanmean(allInvPhaseGainSlopes)
getSEMacrossRows(allInvPhaseGainSlopes)
[p,h,stats]=signrank(allInvPhaseGainSlopes)


saveas(gcf,'bestLogBaseEstimationFig2.png')
%%
fH=figure;subplot(2,1,2)
getJointDistr(1-allFieldDerivWithTimePts,1./allFieldDerivWithTimeVals,linspace(0,1,31),linspace(0,0.008,31),fH);
subplot(2,1,1); getJointDistr(1-allFieldDerivWithTimePts,allFieldDerivWithTimeVals,linspace(0,1,31),linspace(0,1000,31),fH);


%%
figure;
subplot(1,2,1); histogram(allFieldsLogCentralConcOfErr,50)
subplot(1,2,2); histogram(allFieldsLinCentralConcOfErr,50)
close all
numFieldsUsed
setTightSubplots_Spacious

numRhoBins=50;
rhoEdges=linspace(-1,1,numRhoBins+1);
rhoEdges=linspace(0.6,1,numRhoBins+1);
%rhoEdges=linspace(0.8,1,numRhoBins+1);
rhoBins=edgesToBins(rhoEdges);

numErrorBins=50;
numErrorBins=100;
minError=-90;
maxError=90;

minError=-80;
maxError=80;
errorEdges=linspace(minError,maxError,numErrorBins+1);

errorBins=edgesToBins(errorEdges);

figure; 
subplot(2,1,1)
Nlog=histcounts(allFieldLogTimePhaseRs,rhoEdges);
Nlin=histcounts(allFieldTimePhaseRs,rhoEdges);
pLinRho=plot(rhoBins,smooth(Nlin,5)/sum(smooth(Nlin,5)),'k-','LineWidth',4)
hold on
pLogRho=plot(rhoBins,smooth(Nlog,5)/sum(smooth(Nlog,5)),'r-','LineWidth',4)

title(sprintf('R distributions, %d fields',length(allFieldTimeMeanCircErrors) ))
xlim([-1 1])
xlim([0.6 1])
ylim([0 Inf])
legend([pLinRho pLogRho],{'phase-time linear model','phase-log(time) linear model'},'Location','northwest')

xlabel('Correlation coefficient')
%ylabel('Field count')
ylabel('Probability')

box off
subplot(2,1,2)
NlogErr=histcounts(allFieldLogTimeMeanCircErrors,errorEdges);
NlinErr=histcounts(allFieldTimeMeanCircErrors,errorEdges);
pLinErr=plot(errorBins,smooth(NlinErr,5)/sum(smooth(NlinErr,5)),'k-','LineWidth',4)
hold on
pLogErr=plot(errorBins,smooth(NlogErr,5)/sum(smooth(NlogErr,5)),'r-','LineWidth',4)
xlabel('Mean error (deg)')
%ylabel('Field count')
ylabel('Probability')
plot([0 0], ylim,'k--','LineWidth',4)
ylim([0 Inf])


legend([pLinErr pLogErr],{'phase-time linear model','phase-log(time) linear model'},'Location','northwest')
box off

xlim([-30 30])
hold on


title(sprintf('error distributions, %d fields',length(allFieldTimeMeanCircErrors) ))

setFigFontTo(18)
maxFigHalfWidth

 %%%%%%%%%%%%%%%
%STATS
%%%%%%%%%%%%%%%
nanmean(allFieldLogTimePhaseRs)
getSEMacrossRows(allFieldLogTimePhaseRs)
[p,h,stats]=signrank(allFieldLogTimePhaseRs)


saveas(gcf,fullfile(saveFieldPhaseVsLogTimeDir,'linearModelLogTimeFitStatsDistr.png'))
%%

figure; 
subplot(1,2,1)
histogram(allFieldsLogCentralConcOfErr,50)
subplot(1,2,2)
histogram(allFieldsLinCentralConcOfErr,50)

close all
figure;
subplot(1,2,1)
plot(allFieldsLogErrVsTime(:,1),allFieldsLogErrVsTime(:,2),'r.')

subplot(1,2,2)
plot(allFieldsLinErrVsTime(:,1),allFieldsLinErrVsTime(:,2),'k.')

numErrTimeBins=30;

numErrTimeBins=15;
numErrTimeBins=10;
numErrTimeBins=20;
numErrTimeBins=30;
numErrTimeBins=25;
numErrTimeBins=20;
errorTimeEdges=linspace(0,1,numErrTimeBins+1);
errorTimeBins=edgesToBins(errorTimeEdges);

%umErrLevBins=30;
numErrLevBins=90;
errorLevelEdges=linspace(-120, 120,numErrLevBins+1);
figure

%invertedLogTime=invertLogTime(allFieldsLogErrVsTime(:,1),base);
%timeBinsPerRow=discretize(invertedLogTime(:),errorTimeEdges);

%invertedLogTime=invertLogTime(allFieldsLogErrVsTime(:,1),base);
timeBinsPerRow=discretize(allFieldsLogErrVsTime(:,1),errorTimeEdges);

notNaNIdxes=~isnan(timeBinsPerRow) & ~isnan(allFieldsLogErrVsTime(:,2));

%meanLogErrPerTimeBin = splitapply(@nanmean,allFieldsLogErrVsTime(:,2), timeBinsPerRow);
meanLogErrPerTimeBin = splitapply(@circMeanDeg,allFieldsLogErrVsTime(:,2), timeBinsPerRow);

semLogErrPerTimeBin = splitapply(@getSEMacrossRows,allFieldsLogErrVsTime(:,2), timeBinsPerRow);

subplot(1,2,1)
shadedErrorBar(errorTimeBins,meanLogErrPerTimeBin,semLogErrPerTimeBin)

timeBinsPerRow=discretize(allFieldsLinErrVsTime(:,1),errorTimeEdges);
%meanLinErrPerTimeBin = splitapply(@nanmean,allFieldsLinErrVsTime(:,2), timeBinsPerRow);
%semLinErrPerTimeBin = splitapply(@nanstd,allFieldsLinErrVsTime(:,2), timeBinsPerRow);

meanLinErrPerTimeBin = splitapply(@circMeanDeg,allFieldsLinErrVsTime(:,2), timeBinsPerRow);
semLinErrPerTimeBin = splitapply(@getSEMacrossRows,allFieldsLinErrVsTime(:,2), timeBinsPerRow);
%ylim([-40 40])
%figure;
subplot(1,2,2)
shadedErrorBar(errorTimeBins,meanLinErrPerTimeBin,semLinErrPerTimeBin)
%ylim([-40 40])




fH=figure

subplot(1,2,1)
rescaledBehavTime=scaledata(allFieldsLinErrVsTime(:,1),0,0.999);
getJointDistrGivenX(1-rescaledBehavTime,allFieldsLinErrVsTime(:,2),errorTimeEdges,errorLevelEdges,fH)
xlim([0 1])
axis square

hold on
plot(xlim,[0 0],'k--','LineWidth',5)
climits=caxis;

%xlabel('Time to end of field (frac)')

xlabel('Time in field (frac)')
ylabel('Model error (deg)')
title('Phase-time regression error across field')


subplot(1,2,2)
invertedLogTime=timeOffset+invertLogTime( allFieldsLogErrVsTime(:,1),base);
invertedLogTime=scaledata(invertedLogTime,0,0.999);

getJointDistrGivenX(1-invertedLogTime,allFieldsLogErrVsTime(:,2),errorTimeEdges,errorLevelEdges,fH)
xlim([0 1])
axis square
hold on
plot(xlim,[0 0],'k--','LineWidth',5)
%xlabel('Time to end of field (frac)')
xlabel('Time in field (frac)')
ylabel('Model error (deg)')
title('Phase-log(time) regression error across field')
box off

caxis(climits)



setFigFontTo(18)
box off

maxFig
saveas(gcf,fullfile(saveFieldPhaseVsLogTimeDir,sprintf('logTimeVsThetaPhaseErrorOverTimeAllFields.png')))


%xlim([-1 1])

%{
histogram(allFieldLogTimePhaseSlopes,30)
title(sprintf('circSlope distribution (deg/field), %d fields',numFieldsUsed ))

xlim([-360 720])
%}






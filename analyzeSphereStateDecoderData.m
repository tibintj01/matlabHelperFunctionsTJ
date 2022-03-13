
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load data from simulation results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataDir='/Users/tibinjohn/thetaSeq/code/hocFiles/Branco_2010';
close all
nsyn=7;
nsyn=7
directionStr='distToProx';
%directionStr='proxToDist';
%precessionType='log';
precessionType='linear';
allDOF=1

numSims=1;

numStatesPerDim=4

%{
vTraces=readtable(fullfile(dataDir,'tempV.dat'));
for i=0:(nsyn-1)
    gNMDA_Traces.(sprintf('DendSeg%d',i))=readtable(fullfile(dataDir,sprintf('tempNMDA_DendSeg%d.dat',i)));
    gGABA_Traces.(sprintf('DendSeg%d',i))=readtable(fullfile(dataDir,sprintf('tempGABA_DendSeg%d.dat',i)));
end
pspMaxes=readtable(fullfile(dataDir,'tempM.dat'));
%}

numCenters=5;
%numCenters=3;
colorRes=1000;
oversampledColors=jet(colorRes);

colorIDs=linspace(1,colorRes,numCenters);
for i=1:numCenters
    currID=floor(colorIDs(i));
    centerSynapseColors(i,:)=oversampledColors(currID,:);
end

minCenter=ceil(nsyn/2)-floor(numCenters/2)-1;
maxCenter=ceil(nsyn/2)+floor(numCenters/2)-1;

spikeFig3D=figure

if(allDOF)
    minCenter=0;
    maxCenter=0;
end

%for centerSynapseID=1:numCenters
for centerSynapseID=minCenter:maxCenter
    descrString=sprintf('%s_precession_center_neuron_=_%d',precessionType,centerSynapseID);
    
    vTraces=readtable(fullfile(dataDir,sprintf('tempV_CenterSyn%d_%sPrecession.dat',centerSynapseID,precessionType)));
    for i=0:(nsyn-1)
        gNMDA_Traces.(sprintf('DendSeg%d',i))=readtable(fullfile(dataDir,sprintf('tempNMDA_DendSeg%d_CenterSyn%d_%sPrecession.dat',i,centerSynapseID,precessionType)));
        gGABA_Traces.(sprintf('DendSeg%d',i))=readtable(fullfile(dataDir,sprintf('tempGABA_DendSeg%d_CenterSyn%d_%sPrecession.dat',i,centerSynapseID,precessionType)));
    end
    pspMaxes=readtable(fullfile(dataDir,sprintf('tempM_CenterSyn%d_%sPrecession.dat',centerSynapseID,precessionType)));

    dt=0.025; %msec
    saveDt=0.5;
    %simTime=200
        %startTime=50

       simTime=150
    startTime=0
    %timeAxis=dt:dt:simTime;
    timeAxis=(saveDt:saveDt:simTime)-startTime;
    dispTimeLims=[-10 130];
    %dispTimeLims=[-10 150];
    
    
    dispVLims=[-85 -40];
    
    dispVLims=[-85 0];
    %dispVLims=[-75 -69];
    %dispVLims=[-75 -50];
    numTraces=size(vTraces,2)-1;
    maxTotalGgaba=0.1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %set "spike" threshold
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    %thresh=8.17
        %minResponse=7.
    %minResponse=7.9
     minResponse=7.5
    thresh=7.94
     minResponse=0
    thresh=4
    redLineV=-45
    disynapticDelay=0

    statesInSphere=readtable(fullfile(dataDir,'statesInSphere.dat'));
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %detect spike times from inhib conductances
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plotExamples=1;
    plotStep=floor((numStatesPerDim^nsyn)/23);
    plotStep=1;
    if(plotExamples)
        saveCount=0;
        %for simNum=1:plotStep:(numStatesPerDim^nsyn)
        for simNum=1:numSims

            maxgGABA=max(gGABA_Traces.DendSeg0{:,simNum});
            gabaDetectThresh=0.01*maxgGABA;

            inputSpikeTimes=[];

            for i=0:(nsyn-1)
                currGABAtrace=gGABA_Traces.(sprintf('DendSeg%d',i)){:,simNum};
                aboveThreshIdxes=find(currGABAtrace>gabaDetectThresh);
                if(isempty(aboveThreshIdxes))
                   continue 
                end
                firstCrossIdx=aboveThreshIdxes(1);
                firstCrossTime=timeAxis(firstCrossIdx);
                inputSpikeTimes=[inputSpikeTimes; firstCrossTime];
            end
            inputSpikeTimes=inputSpikeTimes-disynapticDelay;
            %figure; 
            %plot(1:length(inputSpikeTimes),inputSpikeTimes,'o')
            %hold on
            %plot([1 length(inputSpikeTimes)],[min(inputSpikeTimes) max(inputSpikeTimes)],'k--')

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %plot data from simulation results
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            figure
            totalGGABA=zeros(size(gGABA_Traces.DendSeg0{:,simNum}));
            axArray=[]
            for i=0:(nsyn-1)
                ax=subplot(nsyn+1,1,i+1)
                axArray=[axArray ax];
                plot(timeAxis,gNMDA_Traces.(sprintf('DendSeg%d',i)){:,simNum},'r','LineWidth',3)
                xlim(dispTimeLims)
                %ylim([0 620])
                   ylim([0 1800])
                totalGGABA=totalGGABA+gGABA_Traces.(sprintf('DendSeg%d',i)){:,simNum};
                ylabel(sprintf('DendSeg%d g_{NMDA}',i))
                hold on
                for j=1:length(inputSpikeTimes)
                    if((i+1)==j)
                        plot([inputSpikeTimes(j) inputSpikeTimes(j)],ylim,'k--','LineWidth',2)
                    end
                end
                box off
            end

            ax=subplot(nsyn+1,1,nsyn+1)
            axArray=[axArray ax];
            plot(timeAxis,vTraces{:,simNum},'k', 'LineWidth',3)
            box off
            ylim(dispVLims)
                xlabel('Time since 1st input (msec)')
                ylabel('Soma V_m (mV)')

                hold on
                plot(xlim,[redLineV redLineV],'r--')
            yyaxis right
            plot(timeAxis,totalGGABA)
               ylim([0 maxTotalGgaba])
               ylabel('total Soma g_{GABA}')
            hold on
            for j=1:length(inputSpikeTimes)
                plot([inputSpikeTimes(j) inputSpikeTimes(j)],ylim,'k--','LineWidth',2)
            end
            maxGGABA=max(totalGGABA);
            %ylim([-maxGGABA maxGGABA])
            xlim(dispTimeLims)
            linkaxes(axArray,'x')
            saveas(gcf,sprintf('dendPosSpikeTimeExcInhibSeqDetection%06d.tif',saveCount))
            saveCount=saveCount+1;
            close all
            %{
            if(simNum==1)
                saveas(gcf,'dendPosSpikeTimeExcInhibSeqDetection.tif')
            elseif(simNum==2)
                saveas(gcf,'dendPosSpikeTimeExcInhibSyncDetection.tif')
            else
                saveas(gcf,'dendPosSpikeTimeExcInhibLogSeqDetection.tif')
            end
            %}
            
        end
    end
    %%
    stateResponses=pspMaxes{:,1};

    numSphereStates=numStatesPerDim^3;
    sphereStateResponses=NaN(numStatesPerDim,numStatesPerDim,numStatesPerDim);
    %count=1;
    %for xi=1:numStatesPerDim
     %   for yi=1:numStatesPerDim
      %      for zi=1:numStatesPerDim

      xDrives=[];
      yDrives=[];
      zDrives=[];
     for i=1:size(statesInSphere,1)     
                xi=round(statesInSphere{i,1}*numStatesPerDim+0.5);
                yi=round(statesInSphere{i,2}*numStatesPerDim+0.5);
                zi=round(statesInSphere{i,3}*numStatesPerDim+0.5);
                 sphereStateResponses(xi,yi,zi)=stateResponses(i);
                 xDrives=[xDrives; xi];
                 yDrives=[yDrives;yi];
                 zDrives=[zDrives;zi];
     end
    %%

    x=1:numStatesPerDim; y=1:numStatesPerDim; z=1:numStatesPerDim;
    %x=unique(x);
    %y=unique(y);
    %z=unique(z);
    plot4D(x,y,z,sphereStateResponses,numStatesPerDim,[minResponse thresh])
    xlabel('(Center neuron - 1) drive')
    ylabel('(Center neuron) drive')
    colormap(jet)
    cb=colorbar
    ylabel(cb,{'Soma PSP max', sprintf('circle radius = (Center neuron + 1) drive (center=%d, edge=%d)',min(zDrives),max(zDrives))})
    %caxis([min(stateResponses) max(stateResponses)])
    caxis([minResponse thresh])
    %caxis([5 10])
    %caxis([4.9 8.2])
    title(sprintf('Dendritic pattern response across input configurations, %s',removeUnderscores(descrString)))
    setFigFontTo(14)

    saveas(gcf,sprintf('allConfigurationsSomaResponse_4D_%s.tif',descrString))
    %%
        
    %thresh=7.5;

    %thresh=7;
    spikingIndices=sphereStateResponses>=thresh;
    spikingMap=zeros(size(sphereStateResponses));
    spikingMap(spikingIndices)=1;
    
    plot4D(x,y,z,spikingMap,numStatesPerDim,[0 1])
    xlabel('(Center neuron - 1) drive')
    ylabel('(Center neuron) drive')
    colormap(jet)
    cb=colorbar
    ylabel(cb,{'isFiring', 'circle radius = (Center neuron + 1) drive'})
    caxis([0 1])
    title({sprintf('Dendritic pattern response across input configurations, %s',removeUnderscores(descrString)),sprintf('Response threshold=%.2f mV',thresh)})
    setFigFontTo(14)
    saveas(gcf,sprintf('allConfigurationsSomaSpikingResponse_%s.tif',descrString))

    firingCoordinates=[];
    %for xi=1:numStatesPerDim
    %    for yi=1:numStatesPerDim
    %        for zi=1:numStatesPerDim
                 for i=1:size(statesInSphere,1)   
                     xi=statesInSphere{i,1}*numStatesPerDim+0.5;
                yi=statesInSphere{i,2}*numStatesPerDim+0.5;
                zi=statesInSphere{i,3}*numStatesPerDim+0.5;
                    if(sphereStateResponses(xi, yi, zi)>thresh)
                        firingCoordinates=[firingCoordinates; xi yi zi];
                    end
                 end

    %        end
    %    end
    %end
    
    figure(spikeFig3D)
    if(isempty(firingCoordinates))
        continue
    end
    
    plot3(firingCoordinates(:,1),firingCoordinates(:,2),firingCoordinates(:,3),'o','Color',centerSynapseColors(centerSynapseID-minCenter+1,:),'MarkerSize',(centerSynapseID+1)*2)
    hold on
    xlim([1 numStatesPerDim])
    ylim([1 numStatesPerDim])
    zlim([1 numStatesPerDim])
    xlabel('(Center neuron - 1) drive')
    ylabel('(Center neuron) drive')
    zlabel('(Center neruon + 1) drive')
    colormap jet
    cb=colorbar
    ylabel(cb,'center neuron for triplet timing modulation')
    caxis([minCenter maxCenter])
    daspect([1 1 1])
    title({sprintf('Dendritic pattern response across input configurations, %s precession',precessionType),sprintf('common response threshold=%.2f mV',thresh)})
    setFigFontTo(14)
    saveas(gcf,sprintf('allConfigurationsSomaSpikingResponse_3D_%s.tif',precessionType))
end
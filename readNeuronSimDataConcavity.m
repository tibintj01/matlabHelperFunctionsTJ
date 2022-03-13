
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load data from simulation results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataDir='/Users/tibinjohn/thetaSeq/code/hocFiles/Branco_2010';
close all
nsyn=7
nsyn=3
directionStr='distToProx';
%directionStr='proxToDist';

vTraces=readtable(fullfile(dataDir,'tempV.dat'));
for i=0:(nsyn-1)
    gNMDA_Traces.(sprintf('DendSeg%d',i))=readtable(fullfile(dataDir,sprintf('tempNMDA_DendSeg%d.dat',i)));
end
pspMaxes=readtable(fullfile(dataDir,'tempM.dat'));

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set parameters for organizing plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt=0.025; %msec
saveDt=0.5;
simTime=200;
startTime=50;
%timeAxis=dt:dt:simTime;
timeAxis=(saveDt:saveDt:simTime)-startTime;

numTraces=size(vTraces,2)-1;

numTrialsPerConcavity=500;
%numConcavitiesTested=12;
numConcavitiesTested=25;
numConcavitiesTested=2;
colormp=jet(numConcavitiesTested)

convexities=linspace(-1,1,numConcavitiesTested);

numSkip=1;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%build soma voltage histograms (time and magnitude) for each sequence convexity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numTs=20;
numTs=100;
%numTs=50;
%tEdges=timeAxis-saveDt; %time bin edges
tEdges=linspace(0,max(timeAxis),numTs); %time bin edges

%numVs=length(timeAxis);
numVs=20;
vMagEdges=linspace(-75,-59,numVs); %voltage bin edges
vMagCenters=edgesToBins(vMagEdges);
tCenters=edgesToBins(tEdges);



%pV_magnitude_time=zeros(length(vMagEdges)-1,length(timeAxis),numConcavitiesTested);
%pV_magnitude_time=zeros(length(vMagEdges)-1,length(tEdges)-1,numConcavitiesTested);
pV_magnitude_time=zeros(length(vMagEdges),length(tEdges),numConcavitiesTested);
for i=1:numTraces
    i
    currConcavityIdx=mod(i-1,numConcavitiesTested)+1;
    if(mod(currConcavityIdx,numSkip)~=0)
        continue
    end
    [N,centers]=hist3([vTraces{:,i} timeAxis(:)], 'Edges',{vMagEdges(:)' tEdges(:)'});
    
    %currQuantizedTrace=discretize(vTraces{:,i},vMagEdges);
   
    %accumlate counts of time-voltage bin pairs
    %for ti=1:length(currQuantizedTrace)
        
        
        %pV_magnitude_time(currQuantizedTrace(ti),ti,currConcavityIdx)= pV_magnitude_time(currQuantizedTrace(ti),ti,currConcavityIdx)+1;
    %end
    pV_magnitude_time(:,:,currConcavityIdx)=pV_magnitude_time(:,:,currConcavityIdx)+N/max(N(:));
end
pV_magnitude_time=pV_magnitude_time/numTrialsPerConcavity;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot gNMDA avg heatmaps (position and time) for each sequence convexity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vTimeHist_fig=figure
for c=1:numConcavitiesTested
    %subplot(7,8,currConcavityIdx)
    c
    subplot(5,5,floor(c/numSkip))
    %omarPcolor(timeAxis,vMagCenters,(pV_magnitude_time(:,:,c)),vTimeHist_fig)
    %omarPcolor(tCenters,vMagCenters,(pV_magnitude_time(:,:,c)),vTimeHist_fig)
     omarPcolor(centers{2},centers{1},(pV_magnitude_time(:,:,c)),vTimeHist_fig)
    
    %imagesc(gNMDA_pos_time(:,:,c))
    xlabel('Time (msec)')
    ylabel('Soma V_m (mV)')
    cb=colorbar
    ylabel(cb,'Probability density')
    caxis([0 1])
    %caxis([0 1250])
    %colormap(parula)
    colormap(jet)
    
     title(sprintf('Seq. convexity=%.2f',convexities(c)))
end
setFigFontTo(14)
maxFig
saveas(gcf,sprintf('voltageTime2Dhistograms_Convexity_%s.tif',directionStr))

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot v traces for each sequence convexity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
for i=1:numTraces
    %skip plots if necessary
    currConcavityIdx=mod(i-1,numConcavitiesTested)+1;
    if(mod(currConcavityIdx,numSkip)~=0)
        continue
    end
    %subplot(7,8,currConcavityIdx)
    subplot(5,5,floor(currConcavityIdx/numSkip))
    plot(timeAxis,vTraces{:,i},'Color',colormp(currConcavityIdx,:),'LineWidth',1)
    hold on
       xlim([-10 150])
     title(sprintf('Seq. convexity=%.2f',convexities(currConcavityIdx)))
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%add colorbars
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
count=0
for i=1:numSkip:numConcavitiesTested
    count=count+1;
    %subplot(7,8,i)
    try
     subplot(5,5,count)
    ylim([-75 -59])
    xlabel('Time (msec)')
    ylabel('Soma V_m (mV)')
    colormap(jet)
    cb=colorbar
    ylabel(cb,'Input seq. convexity')
    caxis([-1 1])
    xlim([-10 150])
    end
end


 maxFig
setFigFontTo(14)
saveas(gcf,sprintf('inputSeqConvexityVsSomaVoltage_%s.tif',directionStr))
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%build gNMDA avg heatmaps (position and time) for each sequence convexity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gNMDA_pos_time=zeros(nsyn,length(timeAxis),numConcavitiesTested);
for c=1:numTraces
    currConcavityIdx=mod(c-1,numConcavitiesTested)+1;
    if(mod(currConcavityIdx,numSkip)~=0)
        continue
    end
    
    for s=1:nsyn
        gNMDA_pos_time(s,:,currConcavityIdx)= gNMDA_pos_time(s,:,currConcavityIdx)+gNMDA_Traces.(sprintf('DendSeg%d',s-1)){:,currConcavityIdx}';
    end
   
end
gNMDA_pos_time=gNMDA_pos_time/numTrialsPerConcavity;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot gNMDA avg heatmaps (position and time) for each sequence convexity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gNMDA_fig=figure
for c=1:numConcavitiesTested
    %subplot(7,8,currConcavityIdx)
    subplot(5,5,floor(c/numSkip))
    omarPcolor(timeAxis,0:(nsyn-1),(gNMDA_pos_time(:,:,c)),gNMDA_fig)
    %imagesc(gNMDA_pos_time(:,:,c))
    xlabel('Time (msec)')
    ylabel('Syn no. (towards soma)')
    cb=colorbar
    ylabel(cb,'gNMDA (pS)')
    caxis([0 1000])
    %caxis([0 1250])
    colormap(parula)
     title(sprintf('Seq. convexity=%.2f',convexities(c)))
end
setFigFontTo(14)
maxFig
saveas(gcf,sprintf('gNMDA_time_pos_heatMaps_Convexity_%s.tif',directionStr))
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot Bezier curve for each sequence convexity to show input shape
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
for i=1:numConcavitiesTested
   	currConcavityFrac=(i-1-numConcavitiesTested/2+0.5)/(numConcavitiesTested/2); %from -1 to 1
   
    %[curveX,curveY] = drawBezier([0 0 currConcavityFrac 1],[1 currConcavityFrac 0 0],100);
    if(currConcavityFrac<0)
        %[curveX,curveY] = drawBezier([0 -currConcavityFrac 1 1],[1 1 -currConcavityFrac 0],100);
           [curveX,curveY] = drawBezier([0 0 1+currConcavityFrac  1],[0 -currConcavityFrac 1 1],100);

    else
        %[curveX,curveY] = drawBezier([0 0 1-currConcavityFrac 1],[1 1-currConcavityFrac 0 0],100);
           [curveX,curveY] = drawBezier([0 currConcavityFrac 1 1],[0 0 1-currConcavityFrac 1],100);

    end
    if(strcmp(directionStr,'proxToDist'))
        plot((1-curveX)*100,curveY*100,'Color',colormp(i,:),'LineWidth',3)
    else
       plot(curveX*100,curveY*100,'Color',colormp(i,:),'LineWidth',3)
    end
    hold on
    ylabel('Relative spike time (msec)')
    xlabel('Distance along dendrite from distal end (um)')
    
end
title('Varying convexity with Bezier curves')
colormap(jet)
    cb=colorbar
    ylabel(cb,'Convexity')
    caxis([-1 1])


 %maxFig
setFigFontTo(18)
saveas(gcf,sprintf('inputSeqConvexityModulationBezierCurves_%s.tif',directionStr))


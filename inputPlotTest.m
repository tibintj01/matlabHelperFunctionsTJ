
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameter settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
%for runID=1:3
for runID=2:3
    plotNormAndDotResponse=0;
    if(runID==1)
          precessionType='linear'
    elseif(runID==2)
        precessionType='log'
    else
        precessionType='normAndDotProductWithWeightSequence';
        plotNormAndDotResponse=1;
    end
    centerSynapseID=0;
   

    dataDir='/Users/tibinjohn/thetaSeq/code/hocFiles/Branco_2010/fluxRawData';
    setTightSubplotsTightLR

    inBatches=1;
    numBatches=21;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %read in simulation result files
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(inBatches)
        statesInSphere=[];
        stateResponses=[];
        for bNum=1:numBatches
            currStatesInSphere=readtable(fullfile(dataDir,sprintf('statesInSphereSection%d.dat',bNum)));
            currStatesInSphere=currStatesInSphere{1:end,:};
            statesInSphere=[statesInSphere; currStatesInSphere];

            if(strcmp(precessionType,'linear'))
              currBatchPspMaxes=readtable(fullfile(dataDir,sprintf('tempM_CenterSyn%d_section%d_%sPrecession.dat',centerSynapseID,bNum,precessionType)));
            elseif(strcmp(precessionType,'log'))
               currBatchPspMaxes=readtable(fullfile(dataDir,sprintf('tempM_CenterSyn%d_section%d_%sPrecession.dat',centerSynapseID,bNum+21,precessionType)));

            end
            currStateResponses=currBatchPspMaxes{1:end,1};
            stateResponses=[stateResponses; currStateResponses(:)];
        end
    else
        statesInSphere=readtable(fullfile(dataDir,'statesInSphere.dat'));
        statesInSphere=statesInSphere{1:end,:};%skip 1st line

        pspMaxes=readtable(fullfile(dataDir,sprintf('tempM_CenterSyn%d_%sPrecession.dat',centerSynapseID,precessionType)));
        stateResponses=pspMaxes{2:end,1};
        size(stateResponses) 
    end

    size(statesInSphere)
    size(stateResponses)
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %get DI of all inputs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    useRepeats=1;
      %useRepeats=0;
    allInputDIs=NaN(size(statesInSphere,1),1);
     inputNoRepeatsFlag=zeros(size(statesInSphere,1),1);
     
     
    for i=1:size(statesInSphere,1)
        currInput=statesInSphere(i,:);
        currInput(currInput<0)=NaN;
        withoutNaNInput=currInput(~isnan(currInput));
        if(length(withoutNaNInput)>1 && length(withoutNaNInput)==length(unique(withoutNaNInput)))
            inputNoRepeatsFlag(i)=1;
        end
        %[~,allInputDIs(i)]=getDI(currInput);
        allInputDIs(i)=getDI(currInput);
    end
    if(useRepeats)
        inputNoRepeatsFlag=ones(size(inputNoRepeatsFlag));
    end

    [~,inputDIsortIDs]=sort(allInputDIs,'descend');
        %figure; plot(allInputDIs(inputDIsortIDs))
        allInputDIs=allInputDIs(inputDIsortIDs);
        statesInSphere=statesInSphere(inputDIsortIDs,:);
        inputNoRepeatsFlag=inputNoRepeatsFlag(inputDIsortIDs);
     if(strcmp(precessionType,'linear') || strcmp(precessionType,'log'))
        stateResponses=stateResponses(inputDIsortIDs);
     end

    %%
    figH=figure
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %plot input configurations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %sentinel for no input
    %statesInSphere(statesInSphere<0)=NaN;
    distFromSoma=1:7;
    inputConfigNum=1:size(statesInSphere,1);

    %testMatrix=randn(size(statesInSphere));

    ax1=subplot(1,2,1)
    %h=pcolor(inputConfigNum,distFromSoma,testMatrix)

    imagesc(fliplr(statesInSphere(:,1:end))')
    %imagesc2(fliplr(statesInSphere(:,1:end))')
    jet2=[1 1 1 ;jet(8)];
    colormap(jet2)
    view([90 90])
    xlabel('Input configuration no.')
    ylabel('CA3-CA1 synapse dendritic distance from soma')
    %set(gca,'yscale','log')
    %cb=colorbar('northoutside')
    cb=colorbar('south')
    leftInc=-1;
    %left bottom width height
    cb.Position=[0.05 0.98 0.47 0.01]

    ylabel(cb,sprintf('normalized CA3 drive (%s theta phase assignment)',precessionType))
    cb.Label.Position=[0.5000 2.3474 0]

    xlim([-Inf Inf])

    %%
    if(plotNormAndDotResponse)
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %plot dot product control
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        inputWeights=[1 2 3 4 5 6 7]';
        inputWeights=inputWeights/norm(inputWeights);
        zeroedStatesInSphere=statesInSphere;
        zeroedStatesInSphere(zeroedStatesInSphere<0)=0;
        normedInputVectors=normMatrixColumns(zeroedStatesInSphere');
        stateResponses=inputWeights'*normedInputVectors;

        maxResponse=37;
        stateResponses=stateResponses*maxResponse;
        
        %%
    end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %plot simulation results
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ax2=subplot(1,2,2)
        plot(stateResponses,'k','LineWidth',0.001)
        hold on
        maxVal=max(stateResponses(:));
        view([90 90])
        xlim([-Inf Inf])
        if(strcmp(precessionType,'log'))
            %ylim([0 58])
            %threshVal=57;  
            %ylim([12 maxVal])
            threshVal=0.975*maxVal;
        else
            %ylim([12 37])
            %threshVal=36;
             %ylim([12 maxVal])
            threshVal=0.975*maxVal;
        end
        %ylim([0 maxVal])
        ylim([12 40])
        plot([1 length(stateResponses)],[threshVal threshVal],'r--','LineWidth',2)
        xticklabels({})

        setFigFontTo(14)

        ylabel('CA1 peak response (mV)')
        linkaxes([ax1 ax2],'x')
        %%
        [peakResponse,peakConfigIdx]=max(stateResponses);
        baseBarHalfWidth=(length(stateResponses)/8)/2;
        baseBarWidth=length(stateResponses)/8;

        isColorTransitionForSynDist=zeros(size(statesInSphere));

        for inputConfigNum=2:size(statesInSphere,1)
            for synDist=1:size(statesInSphere,2)
                if(statesInSphere(inputConfigNum-1,synDist)<statesInSphere(inputConfigNum,synDist))
                    isColorTransitionForSynDist(inputConfigNum,synDist)=1;
                end
            end
        end

        xlimStart=NaN(7,1);
        xlimEnd=NaN(7,1);

        for synDist=1:7
            currStart=peakConfigIdx;
            currEnd=peakConfigIdx;
            while(currStart>0 && isColorTransitionForSynDist(currStart,synDist)==0)
                currStart=currStart-1;
            end
            while(currEnd<size(statesInSphere,1) && isColorTransitionForSynDist(currEnd,synDist)==0)
                currEnd=currEnd+1;
            end
            xlimStart(synDist)=currStart+0.5;
            xlimEnd(synDist)=currEnd-0.5;

        end
        %{

        xlimPrev=1;
        for n=1:5
            currLevelBarWidth=baseBarWidth/(8^n);
            xlimStart(n)=mod(floor(peakConfigIdx/currLevelBarWidth),8)*currLevelBarWidth;

            xlimEnd(n)=xlimStart(n)+currLevelBarWidth;
        end
        %}
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %save all-input simulation results and zooms
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for i=0:7
            if(i==0)
                xlim([-Inf Inf])
            elseif(i<7)
                xlim([xlimStart(i) xlimEnd(i)])
            else
                xlim([xlimStart(7)-0.9 xlimEnd(7)+0.9])
            end
            if(plotNormAndDotResponse)
               title(sprintf('Normalize and dot product, seq. weights (zoom %d)',i))
            else
                title(sprintf('CA1 response to CA3 input configurations (zoom %d)',i))
            end
            saveas(gcf,sprintf('inhibitionGradientModel_ResponseToAllConfigZoom%d_%sPrecession.tif',i,precessionType))
        end
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %save spike triggered input histograms for different thresholds
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %thresholds=[0.975 0.95 0.9 0.875];
         thresholds=[0.975 0.95 0.925 0.90];
        numSlots=7;
        numStates=8;
        spikeTriggeredInputHist=zeros(numStates,numSlots);
        setTightSubplotsTightLR
        setTightSubplotsTight_DI
        hH=figure
        
         diFig=figure
         allInputDIs=allInputDIs(logical(inputNoRepeatsFlag));
         stateResponses=stateResponses(logical(inputNoRepeatsFlag));
         statesInSphere=statesInSphere(logical(inputNoRepeatsFlag),:);
         
        for threshID=1:length(thresholds)
          
            currThresh=thresholds(threshID)*max(stateResponses);
            %spikeIndices=find(stateResponses>currThresh);
            spikeIndices=find(stateResponses>currThresh);
            
            fracInputsThatSpike=length(spikeIndices)/length(stateResponses);
             
            spikeInducingDIs=allInputDIs(spikeIndices);
            figure(diFig)
             subplot(1,length(thresholds),threshID)
            
            %yyaxis left
            nbinsDI=max(allInputDIs);
            %histogram(allInputDIs,nbinsDI)
            [allInputN,edges]=histcounts(allInputDIs,nbinsDI);
            binCenters=edgesToBins(edges);
               plot(binCenters,allInputN/sum(allInputN),'k-','LineWidth',5)
            
            %ylim([0 Inf])
            %yyaxis right
            %if(threshID==length(thresholds))
            %    ylabel('probability')
            %else
            %    %yticklabels({})
            %end
            %histogram(spikeInducingDIs,nbinsDI)
        
            [triggeringInputsN,edges]=histcounts(spikeInducingDIs,edges);
            
         
            hold on
            plot(binCenters,triggeringInputsN/sum(triggeringInputsN),'r-','LineWidth',5)
             title(sprintf('Threshold=%.3f of max response',thresholds(threshID)))
             xlabel('input directionality index')
                  %ylim([0 Inf])
                      legend('all inputs','spike-inducing inputs')
          if(threshID==1)
                ylabel('probability')
                maxY=max([allInputN(:)/sum(allInputN); triggeringInputsN(:)/sum(triggeringInputsN)]);
            else
                yticklabels({})
            end
           % ylim([0 maxY])
           ylim([0 0.425])
             
            for inputI=1:length(spikeIndices)
                currInputID=spikeIndices(inputI);
                currInputPattern=statesInSphere(currInputID,:);
                
                for slotNum=1:numSlots
                    currDrive=max(0,currInputPattern(slotNum));
                    currDriveIdx=ceil(currDrive*numStates+0.5);
                    spikeTriggeredInputHist(currDriveIdx,slotNum)=spikeTriggeredInputHist(currDriveIdx,slotNum)+1;
                end
            end
            spikeTriggeredInputHist=spikeTriggeredInputHist/sum(spikeTriggeredInputHist(:));
            %imagesc(flipud(spikeTriggeredInputHist))
            figure(hH)
              subplot(1,length(thresholds),threshID)
            omarPcolor(1:numSlots,(1:numStates)-1,(spikeTriggeredInputHist),hH)
            if(strcmp(precessionType,'log'))
                %set(gca,'XScale','log')
            end
            colormap(jet)
            if(threshID~=1)
                yticklabels({})
            else
                 ylabel('CA3 neuron drive')
            end
            if(threshID==length(thresholds))
                cbT=colorbar("south")
                ylabel(cbT,'Probability density')
            end
            title({sprintf('Threshold=%.3f of max response',thresholds(threshID)), sprintf('Responds to %.2f%% of all inputs',fracInputsThatSpike*100)})
            xlabel('CA3-CA1 synapse position')
           
               caxis([0 0.1])
            daspect([1 1 1])
            
        end
          %cbT.Label.Position=[0.0500   2.0500         0]
    cbT.Position=[0.3 .875 0.4 0.025]
    
        setFigFontTo(16)
        cbT.Label.FontSize=16
        if(plotNormAndDotResponse)
            uberTitle({sprintf('Robustness of spike-triggering input pattern across thresholds:'), sprintf('%s decoding for all possible inputs',capitalize(precessionType))},20)
        else
            uberTitle({sprintf('Robustness of spike-triggering input pattern across thresholds:'), sprintf('%s delay CA1 decoding for all possible inputs',capitalize(precessionType))},20)
        end
            maxFig

         saveas(hH,sprintf('spikeTriggeringHistogramsAllConfig_%sDecoding.tif',precessionType))
        
         figure(diFig)
             setFigFontTo(16)
        if(plotNormAndDotResponse)
            uberTitle({sprintf('Robustness of spike-triggering input pattern across thresholds:'), sprintf('%s decoding for all possible inputs',capitalize(precessionType))},20)
        else
            uberTitle({sprintf('Robustness of spike-triggering input pattern across thresholds:'), sprintf('%s delay CA1 decoding for all possible inputs',capitalize(precessionType))},20)
        end
            maxFig
            
         saveas(diFig,sprintf('spikeTriggeringDIHistogramsAllConfig_%sDecoding.tif',precessionType))
end
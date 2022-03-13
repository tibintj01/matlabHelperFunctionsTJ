runID=2;


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
        %stateResponseMaxTimes=[];
        wholeTraceResponses=[];
        
        for bNum=1:numBatches
            bNum
            tic
            currStatesInSphere=readtable(fullfile(dataDir,sprintf('statesInSphereSection%d.dat',bNum)));
            currStatesInSphere=currStatesInSphere{1:end,:};
            statesInSphere=[statesInSphere; currStatesInSphere];

            if(strcmp(precessionType,'linear'))
              currBatchPspMaxes=readtable(fullfile(dataDir,sprintf('tempM_CenterSyn%d_section%d_%sPrecession.dat',centerSynapseID,bNum,precessionType)));
            elseif(strcmp(precessionType,'log'))
               currBatchPspMaxes=readtable(fullfile(dataDir,sprintf('tempM_CenterSyn%d_section%d_%sPrecession.dat',centerSynapseID,bNum+21,precessionType)));
               %currBatchWholeTracesTable=readtable(fullfile(dataDir,sprintf('tempV_CenterSyn%d_section%d_%sPrecession.dat',centerSynapseID,bNum+21,precessionType)));
            end
            
            currStateResponses=currBatchPspMaxes{1:end,1};
            currBatchWholeTraces=currBatchWholeTracesTable{:,1:(end-1)};
            stateResponses=[stateResponses; currStateResponses(:)];
            %stateResponseMaxTimes=[stateResponseMaxTimes; currStateResponseMaxTimes(:)];
            %wholeTraceResponses=[wholeTraceResponses; currBatchWholeTraces'];
            toc
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
    
    load(fullfile(dataDir,'logPhasePrecessionAllInputResponseTraces.mat'))

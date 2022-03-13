for chunkNum=1:numChunks
        loadStartCh=1+(chunkNum-1)*loadChunkSize;
        loadEndCh=min(totalNumCh,loadStartCh+loadChunkSize-1);
        
        if(loadStartCh>totalNumCh)
            continue
        end

        disp(sprintf('loading %s ch%d-%d lfp data.....',currSessionName,loadStartCh,loadEndCh))
        
        chIDstr=sprintf('%s_Ch%d',currSessionName,loadEndCh);
            savePrefix=fullfile(saveDir,chIDstr);
            
            if(exist(sprintf('%s_%d-%d_minmax.mat',savePrefix,lowFreq,highFreq),'file'))
                disp('processed files already exist')
                continue
            end
            

        tic
        currLoadChannels=loadStartCh:loadEndCh;
        [Data OrigIndex]=LoadBinary(currDataFile,currLoadChannels);
        toc
        timeAxis=linspace(1/Fs/2,length(Data)/Fs-1/Fs/2,length(Data));

        for ci=1:length(currLoadChannels)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %get cycle information
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            currLFP=Data(ci,:);
            currChNum=currLoadChannels(ci);

            chIDstr=sprintf('%s_Ch%d',currSessionName,currChNum);
            savePrefix=fullfile(saveDir,chIDstr);
            %save filtered LFP
            saveFilteredLFP=0;
             fprintf('saving cycle data for ch %d.....',currChNum)
             tic
            [cycleInfo] = ...
                    getCycleInfoLFP(timeAxis,currLFP, Fs, lowFreq, highFreq, savePrefix,saveFilteredLFP);
                toc
                
                


                 %{
                minTimes=cycleInfo.cycleMinTimes;
                maxTimes=cycleInfo.cycleMaxTimes;

                %thetaAmpZPerCycle=zscoreLFP(cycleInfo.cycleMaxAmp-cycleInfo.cycleMinAmp);
                thetaAmpZPerCycle=(cycleInfo.cycleMaxAmp-cycleInfo.cycleMinAmp);
                %cycleInfo.thetaAmpZPerCycle=thetaAmpZPerCycle;
                
            

                numCycles=length(cycleInfo.cycleMaxTimes);
                [spikeAssignMax spikeAssignMin spikeAssignZero] = assignSpikesToCycles2017(spikeTimes, minTimes, maxTimes);

                mode=2; %phase from max to max
                mode=1;%phase from startMin to endMin
                calcTimeOffsets=1;
                [spikePhase spikeOffsets spikeOffsetsEnd] = assignPhaseToSpikes(spikeTimes, spikeAssignMax, spikeAssignMin, minTimes, maxTimes, mode, calcTimeOffsets, spikeAssignZero);

            allSpikeTimesWithPhases=[spikeTimes(:), spikePhase(:)];
                    %save(
                %}

           
        end %channel loop
    
    end %chunk loop
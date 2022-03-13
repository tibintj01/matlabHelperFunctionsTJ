for chunkNum=1:numChunks
        loadStartCh=1+(chunkNum-1)*loadChunkSize;
        loadEndCh=min(totalNumCh,loadStartCh+loadChunkSize-1);
        
        if(loadStartCh>totalNumCh)
            continue
        end

        disp(sprintf('loading %s ch%d-%d lfp data.....',currSessionName,loadStartCh,loadEndCh))
        
        %chIDstr=sprintf('%s_Ch%d',currSessionName,loadEndCh);
         chIDstr=sprintf('%s_Ch%dRippleInfo.mat',currSessionName,loadEndCh);
            lastInChunkSaveFileName=fullfile(saveDir,chIDstr);
            
            %lastInChunkSaveFileName=sprintf('%s_RippleInfo.mat',savePrefix);
            if(exist(lastInChunkSaveFileName,'file'))
                disp('processed files already exist')
                continue
            end
            

        tic
        currLoadChannels=loadStartCh:loadEndCh;
        [Data OrigIndex]=LoadBinary(currEEGfilePath,currLoadChannels);
        toc
        timeAxis=linspace(1/Fs/2,length(Data)/Fs-1/Fs/2,length(Data));

        for ci=1:length(currLoadChannels)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %get ripple information
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            currLFP=Data(ci,:);
            currChNum=currLoadChannels(ci);

            chIDstr=sprintf('%s_Ch%dRippleInfo.mat',currSessionName,currChNum);
            saveFileName=fullfile(saveDir,chIDstr);
            %save filtered LFP
            saveFilteredLFP=0;
             fprintf('saving ripple data for ch %d.....',currChNum)
             tic
            %[cycleInfo] = ...
            %        getCycleInfoLFP(timeAxis,currLFP, Fs, lowFreq, highFreq, savePrefix,saveFilteredLFP);
             
	     [rippleInfo]=omarRippleDetect(currLFP,Fs);
         save(saveFileName,'rippleInfo')
	     toc
                
           
        end %channel loop
    
    end %chunk loop

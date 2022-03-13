
        
        loadStartCh=min(pyrLayerChannelNums);
        loadEndCh=max(pyrLayerChannelNums);
        
        %{
        if(loadStartCh>totalNumCh)
            continue
        end
        %}

        disp(sprintf('loading %s ch%d-%d lfp data.....',currSessionName,loadStartCh,loadEndCh))
        
        %chIDstr=sprintf('%s_Ch%d',currSessionName,loadEndCh);
            %savePrefix=fullfile(saveLFPcycleInfoDir,chIDstr);
            
            %{
            if(exist(sprintf('%s_fullLFP.mat',savePrefix),'file'))
                disp('processed files already exist')
                continue
            end
            %}

        tic
        currLoadChannels=pyrLayerChannelNums;
        
        [Data OrigIndex]=LoadBinary(currEEGfilePath,currLoadChannels);
        toc
        timeAxis=linspace(1/Fs/2,length(Data)/Fs-1/Fs/2,length(Data));

        rawLFPpyrLayerChs=[];
        rawZscoreLFPpyrLayerChs=[];
        for ci=1:length(currLoadChannels)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %get cycle information
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            currChNum=currLoadChannels(ci);
            rawLFPpyrLayerChs.(sprintf('ch%d',currChNum))=Data(ci,:);
            rawZscoreLFPpyrLayerChs.(sprintf('ch%d',currChNum))=zscoreLFP(Data(ci,:));

            chIDstr=sprintf('%s_Ch%d',currSessionName,currChNum);
            savePrefix=fullfile(saveLFPcycleInfoDir,chIDstr);
          
            disp(sprintf('storing full lfp data for ch %d.....',currChNum))

        end %channel loop
        
        
          save(lfpStructForSessionSavePath,'timeAxis', 'rawLFPpyrLayerChs','rawZscoreLFPpyrLayerChs','pyrLayerChPerShank')
    

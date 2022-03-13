
        
        loadStartCh=min(pyrLayerChannelNums);
        loadEndCh=max(pyrLayerChannelNums);
        
        disp(sprintf('loading %s ch%d-%d lfp data.....',currSessionName,loadStartCh,loadEndCh))
        
        tic
        currLoadChannels=pyrLayerChannelNums;
        
        [Data OrigIndex]=LoadBinary(currEEGfilePath,currLoadChannels);
        toc
        timeAxis=linspace(1/Fs/2,length(Data)/Fs-1/Fs/2,length(Data));

        adjustedCycleInfoPerPyrCh=[];
        
        for ci=1:length(currLoadChannels)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %get cycle information
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            currChNum=currLoadChannels(ci);
            %rawLFPpyrLayerChs.(sprintf('ch%d',currChNum))=Data(ci,:);
            %rawZscoreLFPpyrLayerChs.(sprintf('ch%d',currChNum))=zscoreLFP(Data(ci,:));
            
            disp(sprintf('collecting adjusted cycle lfp data for ch %d.....',currChNum))

            
            currLFP=zscoreLFP(Data(ci,:));
            
            chIDstr=sprintf('%s_Ch%d',currSessionName,currChNum);
            savePrefix=fullfile(saveLFPadjustedCycleInfoDir,chIDstr);
            lowSlowFreq=6;
            highSlowFreq=12;
            
            lowFastFreq=6;
            highFastFreq=40;
            
            saveFilteredLFP=0;
            
            tic
            %{
            [currChAdjustedCycleInfo] = ...
                    getAdjustedCycleInfoLFP(timeAxis,currLFP, Fs, lowSlowFreq, highSlowFreq, lowFastFreq, highFastFreq,savePrefix,saveFilteredLFP);
            %}
                
             [currChAdjustedCycleInfo] = ...
                omarAsymmetryLFP(timeAxis, currLFP, Fs, lowSlowFreq, highSlowFreq, lowFastFreq, highFastFreq, savePrefix)
            toc
                
            
            adjustedCycleInfoPerPyrCh.(sprintf('ch%d',currChNum))=currChAdjustedCycleInfo;

        end %channel loop
        
        
        %save(lfpStructForSessionSavePath,'timeAxis', 'rawLFPpyrLayerChs','rawZscoreLFPpyrLayerChs','pyrLayerChPerShank')
            save(lfpStructForSessionSavePath,'adjustedCycleInfoPerPyrCh','pyrLayerChPerShank')


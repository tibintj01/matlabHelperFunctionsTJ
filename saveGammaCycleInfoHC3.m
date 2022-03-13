  
%lowGammaFreq=30;
%lowGammaFreq=20;
lowGammaFreq=25;
%highGammaFreq=120;
highGammaFreq=150;
        loadStartCh=min(pyrLayerChannelNums);
        loadEndCh=max(pyrLayerChannelNums);

        disp(sprintf('loading %s ch%d-%d lfp data.....',currSessionName,loadStartCh,loadEndCh))
        
        %{
        chIDstr=sprintf('%s_Ch%d',currSessionName,loadEndCh);
            savePrefix=fullfile(saveLFPcycleInfoDir,chIDstr);
            
            if(exist(sprintf('%s_%d-%d_minmax.mat',savePrefix,lowFreq,highFreq),'file'))
                disp('processed files already exist')
                continue
            end
        %}
            

        tic
        currLoadChannels=pyrLayerChannelNums;
        
        [Data OrigIndex]=LoadBinary(currEEGfilePath,currLoadChannels);
        toc
        timeAxis=linspace(1/Fs/2,length(Data)/Fs-1/Fs/2,length(Data));

        for ci=1:length(currLoadChannels)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %get cycle information
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

             
            currLFP=Data(ci,:);
            currChNum=currLoadChannels(ci);

            chIDstr=sprintf('%s_Ch%d',currSessionName,currChNum);
            savePrefix=fullfile(saveLFPcycleInfoDir,chIDstr);
            %save filtered LFP
            saveFilteredLFP=0;
             fprintf('saving cycle data for ch %d.....',currChNum)
             
             try
             tic
            [cycleInfo] = ...
                    getCycleInfoLFP(timeAxis,currLFP, Fs, lowGammaFreq, highGammaFreq, savePrefix,saveFilteredLFP);
                toc
                
             catch ME
                 disp(ME.message)
                 continue
                 
             end
                
           
        end %channel loop
    

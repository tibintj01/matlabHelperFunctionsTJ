function [timeDomainLFPprops] = getFiltAndPhaseLFP(lfp,ch,fpass,saveDir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Description:
%
%
%Input: channel number, lower and upper frequency cutoffs (2 element vector fpass)
% path of directory (appropriately named with subject, session, etc.)
%
%Output:
%
%
%Author: Tibin John, tibintj@umich.edu
%Project directory name: /nfs/turbo/lsa-ojahmed/tibin/segmentedCycleSVDandSpikesProject/scratch 
%Created on 2018-06-24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%extract input variables from struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fs=2000;
chStr=getChStr(ch);

savedPhaseLFPFileName=fullfile(saveDir,sprintf('Ch%sPhaseLFP_lowFreq%.2f_HighFreq%.2f.mat',chStr,fpass(1),fpass(2)));

if(exist(savedPhaseLFPFileName,'file')==2)
        disp('loading phase series.........')
        tic
        phaseLFPData=load(savedPhaseLFPFileName);
        phaseLFP=phaseLFPData.phaseLFP;
        filtLFP=phaseLFPData.filtLFP;
        toc
else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %get filtered lfp
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        filterOrder=2;
        zscoreLFP=1;
        filtLFP=filterLFP(lfp,Fs,fpass(1),fpass(2),filterOrder,zscoreLFP);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %compute lfp phase series 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp('computing phase series.........')
        tic
        hilbertLFP = hilbert(-filtLFP); % NOTE THE INVERSION TO GET 180=MINIMA
        phaseLFP = angle(hilbertLFP);
        % phase ranges from -pi to +pi... make it range from 0 to 360
        phaseLFP = ((phaseLFP/pi) + 1)/2 * 360;
        toc

        disp('saving phase series.........')
        tic
        save(savedPhaseLFPFileName,'phaseLFP','filtLFP')
        toc
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%store output variables in struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timeDomainLFPprops.filtLFP=filtLFP;
timeDomainLFPprops.filtLFP_Z=zscore(filtLFP);
timeDomainLFPprops.phaseLFP=phaseLFP;

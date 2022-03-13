
close all; clear all; clc
dataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles';
showPlots=1;

%sessionNames={'Achilles_10252013','Achilles_11012013','Buddy_06272013', 'Cicero_09102014', 'Cicero_09172014','Gatsby_08282013'};

%Cicero_09172014 has only 1 strongly precessing cell used
sessionNames={'Achilles_10252013','Achilles_11012013','Buddy_06272013', 'Cicero_09102014','Gatsby_08282013'};

halfLineHeight=0.01; %raster line height

for si=1:length(sessionNames)
    currSesName=sessionNames{si};
    trackFilePath=sprintf('%sTrackProps.mat',currSesName);
    currTrackData=load(trackFilePath);
    currTrackData=currTrackData.trackPopInfo;
    %figure;
    rowColInSequenceMatrix=[];
    for di=1:2
        if(strContainsCircSessionName(currSesName))
            if(di==1) %circ only left
                continue
            end
        else
            subplot(2,1,di)
        end
        if(di==1)
           currDirStr='rightward';
        else
           currDirStr='leftward';
        end
        
        if(~isfield(currTrackData,currDirStr))
            continue
        end
        fieldNames=fieldnames(currTrackData.(currDirStr));
        rowColInSequenceMatrix.(currDirStr)=zeros(length(fieldNames),length(fieldNames));
        
        for fi=1:length(fieldNames)
            
            currSpikeTimes=currTrackData.(currDirStr).(fieldNames{fi}).inFieldSpikeTimes;
            currFieldCenterPos=currTrackData.(currDirStr).(fieldNames{fi}).fieldPosCenterM;
            currFieldPosStart=currTrackData.(currDirStr).(fieldNames{fi}).fieldPosStart;
            currFieldPosEnd=currTrackData.(currDirStr).(fieldNames{fi}).fieldPosEnd;
            
             for fii=1:length(fieldNames)
                 otherCellPosStart=currTrackData.(currDirStr).(fieldNames{fii}).fieldPosStart;
                 otherCellPosEnd=currTrackData.(currDirStr).(fieldNames{fii}).fieldPosEnd;
                 
                 hasOverlap=0;
                 if(di==1)
                     %if(otherCellPosStart<currFieldPosEnd || otherCellPosEnd>currFieldPosStart)
                     if(otherCellPosStart<currFieldPosEnd && otherCellPosStart>currFieldPosStart ) 
                         hasOverlap=1;
                     end
                 else
                      %if(otherCellPosStart>currFieldPosEnd || otherCellPosEnd<currFieldPosStart)
                     if(otherCellPosStart>currFieldPosEnd && otherCellPosStart<currFieldPosStart)
                         hasOverlap=1;
                     end
                 end
                 rowColInSequenceMatrix.(currDirStr)(fi,fii)=hasOverlap;
             end
            
            if(showPlots)
                for ti=1:length(currSpikeTimes)
                    %plot([currSpikeTimes(ti) currSpikeTimes(ti)],[currFieldCenterPos-halfLineHeight currFieldCenterPos+halfLineHeight],'k-')
                   if(di==1)
                       plot([currSpikeTimes(ti) currSpikeTimes(ti)],[currFieldPosStart currFieldPosEnd],'k-')
                   else
                       plot([currSpikeTimes(ti) currSpikeTimes(ti)],[currFieldPosEnd currFieldPosStart],'k-')
                   end
                    hold on
                end
            end
            %}
        end
        if(showPlots)
            xlabel('Time (sec)')
            ylabel('Position in track (m)')
            title(currDirStr)
        end
    end
    
    save(trackFilePath,'rowColInSequenceMatrix','currSesName','-append')
    if(showPlots)
        uberTitle(removeUnderscores(currSesName))
        setFigFontTo(18)
        maxFig
    end
end
function []= addFieldResetTimesToDataStruct(filePath)

data=load(filePath);
trackLength=120; %cm
posPerCycleCm=data.spikePerCycleInfo.normPosPerCycle*trackLength;
lapDirectionPerCycle=data.spikePerCycleInfo.lapDirectionPerCycle;
speedPerCycle=data.spikePerCycleInfo.speedPerCycle;

dirFlag=data.spaceTimePhaseInfo.dirFlag;


numCycles=length(data.spikePerCycleInfo.cycleStartTimes);
tElapsedFieldPerCycle=NaN(numCycles,2);
tTotalFieldPerCycle=NaN(numCycles,2);

lastEntryCyclePerCycle=NaN(numCycles,2);
nextLeaveCyclePerCycle=NaN(numCycles,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%***manual start and end only applies to direction manually chosen***%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
entryPositionCm=data.manualFieldStartCm;
leavePositionCm=data.manualFieldEndCm;

middleFieldCm=mean([entryPositionCm leavePositionCm]);

minPosDiff=min(diff(posPerCycleCm));

atMiddleCycleNums=[];
enterFieldCycleNums=[];
leaveFieldCycleNums=[];



%determine which direction of crossing midline corresponds to dirFlag
lapNumPerCycle=data.spikePerCycleInfo.lapNumPerCycle;
middlePositionsThisDir=posPerCycleCm(lapDirectionPerCycle==dirFlag & (lapNumPerCycle==1 | lapNumPerCycle==2));
spatialDiffsThisDir=diff(middlePositionsThisDir);
figure; histogram(spatialDiffsThisDir)
disp('')

figure; plot(posPerCycleCm)
hold on
plot(xlim,[leavePositionCm leavePositionCm],'g')
plot(xlim,[entryPositionCm entryPositionCm],'b')
plot(xlim,[middleFieldCm middleFieldCm],'r')

numLaps=max(lapNumPerCycle);
figure
enterFieldCycleNums=NaN(numLaps,1);
leaveFieldCycleNums=NaN(numLaps,1);

for li=1:numLaps
    currLapCycles=lapNumPerCycle==li;
    if(median(lapDirectionPerCycle(currLapCycles))~=dirFlag)
        %if(median(lapDirectionPerCycle(currLapCycles))==dirFlag)
        continue
     end
  
    
    currLapPositions=posPerCycleCm(currLapCycles);
      currLapCycles(isnan(currLapPositions))=logical(0);
      currLapPositions(isnan(currLapPositions))=[];
      
    firstCycleThisLap=min(find(currLapCycles));
    
    highToLowField=0;
    if((currLapPositions(1)-currLapPositions(end))>0)
        highToLowField=1;
    end
    
    
    if(~highToLowField)
        currSearchCycle=length(currLapPositions);
        while(isnan(currLapPositions(currSearchCycle)) || currLapPositions(currSearchCycle)>entryPositionCm)
            currSearchCycle=currSearchCycle-1;
            if(currSearchCycle<1)
                break
            end
        end
        enterFieldCycleNums(li)=currSearchCycle+firstCycleThisLap-1;
        
        currSearchCycle=1;
        while(isnan(currLapPositions(currSearchCycle)) || currLapPositions(currSearchCycle)<leavePositionCm)
            currSearchCycle=currSearchCycle+1;
            if(currSearchCycle>length(currLapPositions))
                break
            end
        end
        leaveFieldCycleNums(li)=currSearchCycle+firstCycleThisLap-1;
    else
        currSearchCycle=length(currLapPositions);
        while(isnan(currLapPositions(currSearchCycle)) ||currLapPositions(currSearchCycle)<leavePositionCm)
            currSearchCycle=currSearchCycle-1;
            if(currSearchCycle<1)
                break
            end
        end
        enterFieldCycleNums(li)=currSearchCycle+firstCycleThisLap-1;
        
        currSearchCycle=1;
        while(isnan(currLapPositions(currSearchCycle)) || currLapPositions(currSearchCycle)>entryPositionCm)
            currSearchCycle=currSearchCycle+1;
            if(currSearchCycle>length(currLapPositions))
                break
            end
        end
        leaveFieldCycleNums(li)=currSearchCycle+firstCycleThisLap-1;
    end
    %enterFieldCycleNums(li)=
    %leaveFieldCycleNums(li)=
    yyaxis left
    plot(posPerCycleCm(currLapCycles))
    
    
    hold on
    plot(xlim,[leavePositionCm leavePositionCm],'g')
    plot(xlim,[entryPositionCm entryPositionCm],'b')
    try
    p1.Visible='off';
    p2.Visible='off';
    end
    
    p1=plot([leaveFieldCycleNums(li) leaveFieldCycleNums(li)]-firstCycleThisLap+1,ylim,'k','LineWidth',4)
    p2=plot([enterFieldCycleNums(li) enterFieldCycleNums(li)]-firstCycleThisLap+1,ylim,'k--','LineWidth',4)
    
    yyaxis right
    plot(speedPerCycle(currLapCycles))
    
    close all
    
end %lap loop

enterFieldCycleNums(isnan(enterFieldCycleNums))=[];
leaveFieldCycleNums(isnan(leaveFieldCycleNums))=[];
disp('')


for ci=1:numCycles
    if(lapDirectionPerCycle(ci)~=dirFlag)
        continue
    end
    if(posPerCycleCm(ci)<middleFieldCm & posPerCycleCm(ci+1) >=middleFieldCm)
        atMiddleCycleNums=[atMiddleCycleNums; ci];
        currMiddleCycleNum=ci;
        while(posPerCycleCm(currMiddleCycleNum)>entryPositionCm)
            currMiddleCycleNum=currMiddleCycleNum-1;
        end
        enterFieldCycleNums=[enterFieldCycleNums; currMiddleCycleNum];
        
        p1=plot([currMiddleCycleNum currMiddleCycleNum],ylim,'k','LineWidth',4)
        
        currMiddleCycleNum=ci;
        while(posPerCycleCm(currMiddleCycleNum)<leavePositionCm)
            currMiddleCycleNum=currMiddleCycleNum+1;
        end
        leaveFieldCycleNums=[leaveFieldCycleNums; currMiddleCycleNum];
        p2=plot([currMiddleCycleNum currMiddleCycleNum],ylim,'k--','LineWidth',4)
        xlim([0 currMiddleCycleNum])
        disp('')
        
        
    elseif(posPerCycleCm(ci)>middleFieldCm & posPerCycleCm(ci+1) <=middleFieldCm)
        atMiddleCycleNums=[atMiddleCycleNums; ci];
        
        currMiddleCycleNum=ci;
        while(posPerCycleCm(currMiddleCycleNum)>entryPositionCm)
            currMiddleCycleNum=currMiddleCycleNum-1;
        end
        enterFieldCycleNums=[enterFieldCycleNums; currMiddleCycleNum];
        
        currMiddleCycleNum=ci;
        while(posPerCycleCm(currMiddleCycleNum)<leavePositionCm)
            currMiddleCycleNum=currMiddleCycleNum+1;
        end
        leaveFieldCycleNums=[leaveFieldCycleNums; currMiddleCycleNum];
    end
    
    p1.Visible='off';
    p2.Visible='off';
end


for mi=1:length(atMiddleCycleNums)
    currMiddleCycleNum=atMiddleCycleNums(mi);
    
    %while(posPerCycleCm(currCi)>)
        
    %end
    if(lapDirectionPerCycle(ci)<0)
        currDI=1;
    else
        currDI=2;
    end
    
    currPos=posPerCycleCm(ci);
    
    
    %lastEntryCyclePerCycle(ci,currDI)=
end


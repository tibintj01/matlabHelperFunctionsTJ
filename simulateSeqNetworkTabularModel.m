close all
if(~exist('statesInSphere','var'))
    clear all; clc
    loadTabularCA1Model
end

size(statesInSphere)
size(stateResponses)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set up random feedforward network
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nCA3=7;
nCA1=5040;
postPreConn=zeros(nCA1,nCA3);

Vr=-75;

nCA1_DIs=NaN(nCA1,1);

%each CA3 neuron goes to around 1000 CA1 cells (at approximately the same
%dendritic layer

for post_i=1:nCA1
    %posIsTaken=zeros(nCA3,1);
    postPreConn(post_i,:)=randsample(7,nCA3);%random number from 1 to 7 no repeats
    %for pre_i=1:nCA3
        %postPreConn(post_i,pre_i)=randi(8,1)-1; %random number from 0 to 7
        %postPreConn(post_i,randi(nCA3,1))=pre_i; %random number from 0 to 7
        
    %end
    %nCA1_DIs(post_i)=getDI(postPreConn(post_i,:));
end

%{
divNum=3000;
for pre_i=1:nCA3
    post_Is=randi(nCA1,divNum,1);
    postPreConn(post_Is,pre_i)=pre_i+randi(3,1)-2; %random number from 0 to 7
end
postPreConn(postPreConn<0)=0;
postPreConn(postPreConn>7)=7;

%}
    
for post_i=1:nCA1
    nCA1_DIs(post_i)=getDI(postPreConn(post_i,:));
end
[~,sortI]=sort(nCA1_DIs);

postPreConn=postPreConn(sortI,:);
figure
subplot(2,1,1)
imagesc(postPreConn')
subplot(2,1,2)
histogram(nCA1_DIs)
xlim([1 28])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%simulate CA1 recruitment stats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CA3 sequence simulated 

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%simulate activity bump movement response
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sentinelVal=-0.1250;
numCycles=12;
numCycles=14;

activitySeq=zeros(numCycles+1,nCA3);

activityBump=1:7;

cycleCount=1;
for shift=-numCycles/2:numCycles/2
    shiftedBump=circshift(activityBump,shift)
    startIdx=max(1,1+shift);
    endIdx=min(nCA3,nCA3+shift);
    
    activitySeq(cycleCount,startIdx:endIdx)=shiftedBump(startIdx:endIdx);
    cycleCount=cycleCount+1;
end

%activitySeq=(activitySeq+0.5)/8;
%activitySeq(activitySeq<0.1)=sentinelVal;
%activitySeq

%posMatchResponse=NaN(numCycles,1);
patternSensitiveResponse=[];
cycleLength=150;
cycleLength=100;
cycleLength=125;
peakTimeFracPerCycle=NaN(numCycles+1,1);
amplificationThresh=-50;
for cycle=1:(numCycles+1);
    currActivityPattern=activitySeq(cycle,:);
    
    currVtrace=sequenceToTableFullResponse(currActivityPattern,statesInSphere,wholeTraceResponses);
    [peakVal,maxIdx]=max(currVtrace);
    if(peakVal>amplificationThresh)
        peakTimeFracPerCycle(cycle)=maxIdx/cycleLength;
    end
    patternSensitiveResponse=[patternSensitiveResponse; currVtrace(:)];
    %dotProds=statesInSphere*currActivityPattern(:);
    %[~,bestMatchIdx]=max(dotProds);
    %posMatchResponse(cycle)=stateResponses(bestMatchIdx);
end
timeAxis=(1:length(patternSensitiveResponse))/2;
figure
%plot(currVtrace)
plot(timeAxis,patternSensitiveResponse)
hold on
for cycleEnd=cycleLength/2:cycleLength/2:length(patternSensitiveResponse)/2
    plot([cycleEnd cycleEnd],ylim,'k--','LineWidth',3)
end
plot(xlim,[amplificationThresh amplificationThresh],'r--','LineWidth',5)
ylabel('V (mV)')
xlabel('Time (msec)')
title('CA3 bump movement plus CA1 sequence sensitivity induced phase precession')
setFigFontTo(18)
saveas(gcf,'activityBumpMovementInducedSubthresholdResponses.tif')
figure
plot(peakTimeFracPerCycle,'bo','MarkerSize',6)
hold on
plot(peakTimeFracPerCycle,'-','LineWidth',4)
ylim([0 1])
xlabel('cycle')
ylabel('CA1 subthreshold peak time (fraction of theta cycle)')
title('CA3 bump movement plus CA1 sequence sensitivity induced phase precession')
setFigFontTo(18)
saveas(gcf,'activityBumpMovementInducedPhasePrecession.tif')

figure; omarPcolor(1:size(activitySeq,2), 1:size(activitySeq,1),activitySeq)
ylabel('Theta cycle number')
xlabel('CA1 dendritic position')
cb=colorbar
ylabel(cb,'CA3 input theta phase')
setFigFontTo(18)
title('CA3 activity bump movement')
saveas(gcf,'activityBumpMovement.tif')
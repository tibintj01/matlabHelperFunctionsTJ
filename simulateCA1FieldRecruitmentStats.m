close all
setTightSubplotsTight_Recruitment
if(~exist('statesInSphere','var'))
    clear all; clc
    loadTabularCA1Model
end

for r=1:2
    if(r==1)
        randNetworkBased=1;
    else
        randNetworkBased=0;
    end

    if(randNetworkBased)
        networkStr='Random';
    else
       networkStr='Structured';
    end

    rngSeed=4;
    rngSeed=2;
    rngSeed=1;
    rngSeed=3;

    rng(rngSeed)

    size(statesInSphere)
    size(stateResponses)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %set up random feedforward network
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nCA3=10;
    nCA3=15;
    %nCA3=20;
    %nCA3=25;
    %nCA3=50;
    nCA1=5000;
    nCA1=500;
    nCA1=1000;
    nCA1=1500;
    nCA1=1250;
    nCA1=1000;
    nCA1=500;
    nCA1=250;
    nCA1=200;
    nCA1=150;
    nCA1=1000;
    nCA1=500;
    nCA1=300;
    nCA1=250;
    nCA1=1000;
    nCA1=750;
    nCA1=500;
    nCA1=300;
    nCA1=250;
    postPreConn=zeros(nCA1,nCA3);

    Vr=-75;

    nCA1_DIs=NaN(nCA1,1);

    %each CA3 neuron goes to around 1000 CA1 cells (at approximately the same
    %dendritic layer

        for post_i=1:nCA1
            %posIsTaken=zeros(nCA3,1);
            %postPreConn(post_i,:)=randsample(7,nCA3,1);%random number from 1 to 7 no repeats
            % if(seqBased)
                 rand1=randsample(nCA3,7,1);
                 rand2=randsample(7,7,0);
             %else
              %  rand1=randsample(nCA3,7,1);
               % rand2=repelem(3.5,7);
             %end
            postPreConn(post_i,rand1)=rand2;%random number from 1 to 7 no repeats

        end
    if(~randNetworkBased)
        connMean=mean(postPreConn(:));
        %connStd=std(postPreConn(:));
        %postPreConn=normrnd(connMean,connStd,[nCA1 nCA3]);
        %postPreConn=postPreConn(randperm(nCA1),randperm(nCA3));
        structuredGroupSize=floor(nCA1/nCA3)
        postPreConn=zeros(size(postPreConn));
        for groupNum=1:nCA3
            startIdx=(groupNum-1)*structuredGroupSize+1;
            endIdx=startIdx+structuredGroupSize;
            postPreConn(startIdx:endIdx,groupNum)=normrnd(5,1,length(startIdx:endIdx),1);

        end

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
        currSeq=postPreConn(post_i,:);
        nCA1_DIs(post_i)=getDI(currSeq(currSeq>0));
    end
    [~,sortI]=sort(nCA1_DIs);

    %postPreConn=postPreConn(sortI,:);
    %%
    fH=figure
    subplot(2,3,4)
    imagesc(postPreConn')
    colormap(gca,jet)
    colorbar
    title('Connectivity')
    xlabel('CA1 cell no.')
    ylabel('CA3 cell no.')
    %subplot(2,3,3)
    %histogram(nCA1_DIs)
    %xlim([1 28])

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %simulate CA1 recruitment stats
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %CA3 sequence simulated as random selection of 7 cells
    %updated by increasing sequence by one (phase precession)
    %and replacing maxed-out cell by one random CA3 cell
    %and repeating ad infinitum; 
    %-->recruitment stats? Field number distribution?



    ca3ActivitySeries=getRandomCA3Trajectory(nCA3,rngSeed,randNetworkBased);

    [ca1ActivationSequence,strengths]=getCA1Sequence(ca3ActivitySeries,postPreConn);


    occurenceCounts=zeros(nCA1,1);

    numDistSteps=length(ca1ActivationSequence);
    recruitmentFracSoFar=zeros(numDistSteps,1);
    fieldRasterMatrix=NaN(numDistSteps,nCA1);
    fieldRasterMatrix=zeros(numDistSteps,nCA1);
    for i=1:numDistSteps
        occurenceCounts(ca1ActivationSequence(i))=occurenceCounts(ca1ActivationSequence(i))+1;
        recruitmentFracSoFar(i)=sum(occurenceCounts>0)/nCA1;
        fieldRasterMatrix(i,ca1ActivationSequence(i))=1;
    end


    [counts,sortIdx]=sort(occurenceCounts);
    fieldRasterMatrix=fieldRasterMatrix(:,sortIdx);

    %[counts,sortIdx]=sort(ca1ActivationSequence(occurenceCounts));

    subplot(2,3,[1 2 3])
    omarPcolor((1:size(fieldRasterMatrix,1))/size(fieldRasterMatrix,1),1:size(fieldRasterMatrix,2),fieldRasterMatrix',fH)
    colormap(gca,copper)
    cb=colorbar
    xlabel('Position (a.u.)')
    ylabel('CA1 cell no.')
    ylabel(cb,'Place field activity')

    title(sprintf('CA1 place cell recruitment vs distance: %s CA3 feedforward network',networkStr))

    %cellNums=1:nCA1;
    %cellRemapping=sortIdx(cellNums);
    %plot(sortIdx(ca1ActivationSequence),'ko','MarkerSize',3)

    %figure
    %subplot(2,1,1)
    %histogram(occurenceCounts,0:35)
    subplot(2,3,5)
    histogram(occurenceCounts,0:(max(occurenceCounts)+1))
    xlabel('Fields per CA1 place cell')
    ylabel('Number of CA1 cells')

    %figure;
    subplot(2,3,6)
    plot((1:length(recruitmentFracSoFar))/length(recruitmentFracSoFar),recruitmentFracSoFar)
    ylim([0 1])
    xlabel('Position (a.u.)')
    ylabel('Proportion of cells recruited')
    %%

    axes('Position',[0.875 .385 .1 .1])
    box on
    %subplot(2,3,4)
    plot(recruitmentFracSoFar)
    ylim([0 1])
    set(gca,'xscale','log')

    setFigFontTo(18)
    maxFig

    saveas(gcf,sprintf('%sNetworkStructure_CA1recruitmentVsDistanceStats.tif',networkStr))
end


function [] = checkSpeed(currFilePath)
%check is speed fields make sense
maxSpeed=1.5; %m/s - artifact?

data=load(currFilePath);

speedFieldNames={'oneDspeedPerTimeStep','absSpeedPerTimeStep','filledSignedSpeedPerTime','speedPerTimeStep'}
   
%figure

for i=1:length(speedFieldNames)
    %subplot(2,2,i)
    %plot(data.positionTimeAxis,  data.(speedFieldNames{i}))
    
    artifactIdxes=(data.(speedFieldNames{i})>maxSpeed | data.(speedFieldNames{i})<-maxSpeed );
    data.(speedFieldNames{i})(artifactIdxes)=NaN;
    
end

oneDspeedPerTimeStep=data.oneDspeedPerTimeStep;
absSpeedPerTimeStep=data.absSpeedPerTimeStep;
filledSignedSpeedPerTime=data.filledSignedSpeedPerTime;
speedPerTimeStep=data.speedPerTimeStep;

save(currFilePath,'oneDspeedPerTimeStep','absSpeedPerTimeStep','filledSignedSpeedPerTime','speedPerTimeStep','-append')

%data=load(currFilePath);


%{
figure

for i=1:length(speedFieldNames)
    subplot(2,2,i)
    plot(data.positionTimeAxis,  data.(speedFieldNames{i}))
    
    artifactIdxes=(data.(speedFieldNames{i})>maxSpeed);
    data.(speedFieldNames{i})(artifactIdxes)=NaN;
    
end
%}
disp('')


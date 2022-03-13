close all; clear all; clc
%%%%%%%%%%%%%%%%%%%%%
%test case #1:
%linear relationship 
%with gaussian noise
%%%%%%%%%%%%%%%%%%%%%
distNorm=0:0.001:1';
%noise=50*(rand(length(distNorm),1)-0.5);

for ni=1:2
    
    if(ni==1)
        gaussStdev=10;
        noiseStr='low';
    elseif(ni==2)
        gaussStdev=45;
        noiseStr='high';
    end


    noise=gaussStdev*(normrnd(0,1,length(distNorm),1));

    shapeNames={'linear','sublinear','supralinear'};

    for si=1:length(shapeNames)
        if(si==1)
            phasesDeg=(1-(distNorm(:)))*360;
        elseif(si==2)
            [~, convexIncreasing] = getBezierGivenConvexity(0.8,length(distNorm));
            phasesDeg=fliplr(convexIncreasing*360);
        elseif(si==3)
            [~, concaveIncreasing] = getBezierGivenConvexity(-0.8,length(distNorm));
            phasesDeg=fliplr(concaveIncreasing*360);
        end

        phasesDeg=phasesDeg(:)+noise(:);
        phasesDeg=mod(phasesDeg,360);

        distBinWidth=0.05;

        titleStr= sprintf('%s_%s_noise',shapeNames{si},noiseStr) ;
        [devFromLinearVsSpace] = ...
            getLinearDeviationOfCylDataTrueOrtho(phasesDeg,distNorm,distBinWidth,titleStr);
    end
    
end

%figure; plot(devFromLinearVsSpace)
%phasesDeg=72*(log(1-distNorm(:))+5)+noise(:);

%phasesDeg=(1-(distNorm(:)).^8)*360+noise(:);
%phasesDeg=(1-(distNorm(:)).^8)*360+noise(:).^8;
%phasesDeg=(1-(distNorm(:)))*300+noise(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%test case #2:
%superlinear relationship 
%with gaussian noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%test case #3:
%sublinear relationship 
%with gaussian noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%
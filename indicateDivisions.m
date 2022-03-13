function [] = indicateDivisions(divBoundaries,includeText)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Description:
%
%
%Input:
%
%
%Output:
%
%
%Author: Tibin John, tibintj@umich.edu
%Project directory name: /nfs/turbo/lsa-ojahmed/tibin/spikeDynamicsAnalysisTibin/checkSeizureSortingQuality 
%Created on 2019-01-06
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(~exist('includeText'))
	includeText=1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%extract input variables from struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on

ylims=ylim;
xlims=xlim;

minY=ylims(1);
maxY=ylims(2);

minX=xlims(1);
maxX=xlims(2);

divBoundaries=divBoundaries-divBoundaries(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compute output variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('computing indicateDivisions fxn output.........')

%divColors={'k','m'}
%divColors={[200 200 200]/255,'m'}
divColors={[200 200 200]/255,getGrayRGB()}

if(includeText)
	%plot([minX 0],[maxY maxY],'-','Color',divColors{1},'LineWidth',8)
	textX=mean([minX 0]);
	textY=maxY*0.97;
	text(textX,textY,sprintf('%s','pre-seizure'))

	textX=mean([divBoundaries(end) maxX]);
	textY=maxY*0.97;
	text(textX,textY,sprintf('%s','post-seizure'))

	textX=minX;
	textY=maxY*0.97;
	text(textX,textY,sprintf('%s','SZ Div:'))
end

for i=1:length(divBoundaries)
        %plot([divBoundaries(i) divBoundaries(i)],ylim,'k--','LineWidth',2)
	xtra=(maxY-minY)*0;
        plot([divBoundaries(i) divBoundaries(i)],[minY maxY+xtra],'--','Color',getGrayRGB(),'LineWidth',2)

        colorIdx=mod(i,2)+1;

        if(i<length(divBoundaries))
                %plot([divBoundaries(i) divBoundaries(i+1)],[maxY maxY],'-','Color',divColors{colorIdx},'LineWidth',8)
		if(includeText)
			textX=mean([divBoundaries(i) divBoundaries(i+1)]);
			textY=maxY*0.97;
			text(textX,textY,sprintf('%d',i))
        	end
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%store output variables in struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = plotVertLine(xPos,colLetter)
	hold on
	plot([xPos xPos],ylim,sprintf('%s--',colLetter),'LineWidth',2.5)

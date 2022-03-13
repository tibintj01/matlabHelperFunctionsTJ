function [] = addColorBarRight(cbLabel)

	axPos=get(gca,'Position');

	leftEdge=axPos(1)+axPos(3);
	bottomEdge=axPos(2)+0.02;
	

	width=0.01;
	height=axPos(4);
	cb=colorbar;
	set(cb,'Units','normalized','position',[leftEdge bottomEdge width height]);
	ylabel(cb,cbLabel)


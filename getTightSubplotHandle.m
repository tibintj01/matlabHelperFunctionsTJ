function [tightSubplotH]=getTightSubplotHandle(numRowsSubP,numColsSubP,plotGapWidth,plotGapHeight,plotWidthMargins,plotHeightMargins)
%  Example (from tight_subplot.m): ha = tight_subplot(3,2,[.01 .03],[.1 .01],[.01 .01])
%           for ii = 1:6; axes(ha(ii)); plot(randn(10,ii)); end
%           set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')

if(~exist('plotGapWidth'))
	plotGapWidth=0.01;
end
if(~exist('plotGapHeight'))
	plotGapHeight=0.01;
end

%plotHeightMargins=[0.1 0.01];
%plotHeightMargins=[0.01 0.1];

if(~exist('plotHeightMargins'))
	plotHeightMargins=[0.01 0.11];
end
if(~exist('plotWidthMargins'))
	plotWidthMargins=[0.01 0.01];
end
%plotHeightMargins=[0.15 0.15];
%plotWidthMargins=[0.1 0.1];

tightSubplotH=tight_subplot(numRowsSubP,numColsSubP,[plotGapHeight plotGapWidth],plotHeightMargins,plotWidthMargins)



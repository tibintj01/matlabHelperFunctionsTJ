function [cellIdxWithinChan]=getCellIdxWithinChannel(totalCellIdx,cellChannels)

	cellChannel=cellChannels(totalCellIdx);
	firstChanOccurence=getIdxOf(cellChannel,cellChannels);

	cellIdxWithinChan=totalCellIdx-firstChanOccurence+1;

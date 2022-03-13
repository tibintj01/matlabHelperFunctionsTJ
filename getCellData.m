function [spikeTimes,waveforms,chanNum]=getCellData(ptID,cellNum)
 		chanNum=cellChannels(cellNum);
                nexFileName=getNexFileName(chanNum);
                cellIdxInNex=getCellIdxWithinChannel(cellNum,cellChannels);
                nexData=readNexFile(nexFileName);

                cellSortedSpikes(cellNum).spikeTimes=nexData.neurons{cellIdxInNex}.timestamps;
                %each row is spike, columns represent increasing time
                cellSortedSpikes(cellNum).waveForms=nexData.waves{cellIdxInNex}.waveforms';

                cellSortedSpikes(cellNum).chanNum=chanNum;


function [] = plotOriginalSpikeTimes(unitData)
positionTimeAxis=unitData.positionTimeAxis;

originalSpikeTimes=unitData.unitInfo.spikeTimes;

                        originalSpikeTimes=originalSpikeTimes(originalSpikeTimes>=min(positionTimeAxis)...
                                                        & originalSpikeTimes<=max(positionTimeAxis));
                        plot(  originalSpikeTimes, repmat([1],size( originalSpikeTimes)),'r*','MarkerSize',10) 
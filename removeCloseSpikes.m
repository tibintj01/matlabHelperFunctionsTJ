function [sparseSpikeTimes]=removeCloseSpikes(spikeTimes,minISI)
                idxNextSpikeTooClose=[];
                ISIs=[diff(spikeTimes) Inf];
                for i=1:length(ISIs)
                        if(ISIs(i)<=minISI)
                                idxNextSpikeTooClose=[idxNextSpikeTooClose i];
                        end
         
                end

                spikeTimes(idxNextSpikeTooClose)=[];
                sparseSpikeTimes=spikeTimes;

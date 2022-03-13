
                figure(1)
                waveforms=data.avgwave{cellNum}';
                subplot(numPanels,numSessions,1:numSessions)
                plot(waveforms(:))
                %plot(Spikes{cellNum},cellNum*ones(length(Spikes{cellNum}),1),'o' )
                %figure; plot( ThPrecess{cellNum,1}(:,1),ThPrecess{cellNum,1}(:,2),'o')

                title({sprintf('Cell %d',cellNum),removeUnderscores(fileNameRoot)})



                %for sessionID=1:numSessions
                    %subplot(numPanels,numSessions,numSessions+sessionID)
                    subplot(numPanels,numSessions,numSessions+1)
                    if(size(data.ratemap{cellNum,sessionID},1)==size(data.ratemap{cellNum,sessionID},2))
                        rates=data.ratemap{cellNum,sessionID};
                        imagesc(rates)
                        title(sprintf('max rate: %.1f Hz',max(rates(:))))
                    else
                        posAxis=linspace(0,maxDist,length(data.ratemap{cellNum,sessionID}));
                        plot(posAxis,data.ratemap{cellNum,sessionID})
                        xlim([0 maxDist])
                        xlabel('Distance (cm)')
                        ylabel('Firing rate (Hz)')
                    end
                    hold on
                %end
                %for sessionID=1:numSessions
                %for sessionID=1:1

                   try
                        phases=data.ThPrecess{cellNum,sessionID}(:,2);
                 
                    subplot(numPanels,numSessions,numSessions*2+1)
                    position=data.ThPrecess{cellNum,sessionID}(:,1)*maxDist;    
                    plot(position,phases,'ko','MarkerSize',1)
                    hold on
                    plot(position,phases+360,'ko','MarkerSize',1)
                    xlabel('Distance (cm)')
                    ylabel('Theta phase (deg)')
                    xlim([0 maxDist])
                    ylim([0 720])
                    
                   catch ME
                       disp(ME)
                       close all
                       %continue
                   end
                    phSlope=data.measures(cellNum,14,sessionID);
                    phRSquare=data.measures(cellNum,15,sessionID);
                    title({sprintf('phPrecess slope: %.1f',phSlope),sprintf('R^2=%.3f',phRSquare)})
                    hold on
                %end
                    drawnow
                    maxFig
                    setFigFontTo(16)
                    saveas(gcf,fullfile(saveDir,sprintf('%s_cell%d_ses%d.tif',fileNameRoot,cellNum,sessionID)))
                    close all
if(di==1)
	plot(flipud(spaceTimePhaseInfo.placeField(:)),'k-')
	fieldStart=trackLength-spaceTimePhaseInfo.startFieldCm;
	fieldEnd=trackLength-spaceTimePhaseInfo.endFieldCm;
else
	plot(spaceTimePhaseInfo.placeField(:),'k-')
	fieldStart=spaceTimePhaseInfo.startFieldCm;
	fieldEnd=spaceTimePhaseInfo.endFieldCm;
end
xlim([0 trackLength])
xlabel('Distance along track (cm)')
ylabel('Firing rate (Hz)')
hold on
plot([currStartCm currStartCm],ylim,'b','LineWidth',1)
plot([currEndCm currEndCm],ylim,'r','LineWidth',1)
title('Spike rate vs Distance')

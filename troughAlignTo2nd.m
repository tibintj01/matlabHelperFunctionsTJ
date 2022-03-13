function [aligned] = troughAlignTo2nd(wav1,wav2)

	[dummy,wav1bin]=min(wav1);
	[dummy,wav2bin]=min(wav2);

	aligned=shiftVals(wav1,wav2bin-wav1bin);
	aligned=aligned(:);

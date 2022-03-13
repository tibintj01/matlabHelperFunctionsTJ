function [] = detectCessationTime(firingRateOverTime,timeAxis)


	%return bin and time corresponding to latest bin such that all bins afterwards
	%last bin at 20% of max firing rate - check each visually	

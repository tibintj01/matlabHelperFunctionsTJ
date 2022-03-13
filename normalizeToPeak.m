function [peakNormed]=normalizeToPeak(values)
	peakNormed=values/max(values);

function [wave]=getSineWave(t,freq,amp)
	wave = amp * sin(2 * pi * freq * t);


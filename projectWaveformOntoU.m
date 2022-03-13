function [pc1,pc2]=projectWaveformOntoU(currWaveform,U)
	%U is a matrix with columns the singular vectors of a group waveforms
	%pc1=dot(currWaveform(:),U(:,1))/sigma(1);
	%pc2=dot(currWaveform(:),U(:,2))/sigma(2);
	pc1=dot(currWaveform(:),U(:,1));
	pc2=dot(currWaveform(:),U(:,2));
	

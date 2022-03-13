% Spk = LoadSpk(FileName, nChannels, SpkSampls, SpikesToLoad)
%
% Loads a .spk file
%
% returns a 3d array Spk(Channel, Sample, Spike Number)
%
% nChannels will default to 4
% SpkSampls will default to 32
% Spikes2Load will default to [] - i.e. load them all,except noise

function [Spk,nChannels] = LoadSpk(FileName, nChannels, SpkSampls, Spikes2Load)


if (nargin<2) nChannels = 8; end;
if (nargin<3) SpkSampls = 32; end;
if (nargin<4) Spikes2Load = inf; end;

correctNumChannels=0;
while(correctNumChannels==0)
    Spk = bload(FileName, [nChannels, Spikes2Load*SpkSampls],0,'short=>single');

    nSpikes = size(Spk, 2)/SpkSampls;
    %%Tibin Jan 11, 2020
    if(ceil(nSpikes)~=floor(nSpikes))
        %error('wrong number of channels!!!') %check xml for correct number of channels in this shank (8, 10, etc.)
        disp('wrong number of channels!!!') %check xml for correct number of channels in this shank (8, 10, etc.)
        nChannels=nChannels-1;
        disp(sprintf('trying %d channels',nChannels))
    else
        correctNumChannels=1;
    end
end

Spk = reshape(Spk, [nChannels, SpkSampls, nSpikes]);


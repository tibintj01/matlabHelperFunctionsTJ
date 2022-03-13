close all ; clear all; clc

t=0:0.001:10; %sec

asymTheta=getAsymWave(2*pi*t);

%%%%%%%%%%%%%%%%%%%%%
%plots
%%%%%%%%%%%%%%%%%%%%%
figure;
plot(t,asymTheta,'k-')
t=linspace(0,5,1000);
f0=7;
%f=1.1;
f=linspace(7,7.5,1000);

refOsc=sin(2*pi*f0*t);
modulatedOsc=sin(2*pi*f.*t);

figure

subplot(2,1,1)
plot(t,refOsc)
hold on
plot(t,modulatedOsc+2)
title({'Speed as a squeeze mapping of internal time and space',sprintf('f=%.2f',f(end)),'time dilation'})
xlabel('Time (sec)')
subplot(2,1,2)

[upEnv,lowEnv]=envelope(refOsc+modulatedOsc);
plot(t,upEnv)
hold on
plot(t,lowEnv)
title('length contraction')
disp('')
xlabel('Length')


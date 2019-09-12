clc;
clear all;
close all;

M=100;
h=2*pi/M;
x=-pi:h:pi;
y=sawtooth(x);
y = sin(2*x)+sin(3*x)+sin(4*x) + 6*cos(6*x)+ 7*cos(7*x);
L=length(y);

f = 2*fft(y(1:L-1))/(L-1);
a=real((f));
b=imag((f));

b(1:2:end) = -b(1:2:end);

figure(1)
plot(abs(f),'-*')

figure(2)
plot(x,y)



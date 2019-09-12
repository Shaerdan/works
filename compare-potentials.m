clc;
clear all;
close all;

x = -1:0.01:1
r = 0.5;
x1 = 0.7;
x2 = 0.75;
y1 = log(1+r.^2 -2*r*cos(x-x1)) - log(1+r.^2 -2*r*cos(x-x2));
y2 = 1./(1+r.^2 -2*r*cos(x-x1)) - 1./(1+r.^2 -2*r*cos(x-x2));

figure(1)
plot(x,y1);

figure(2)
plot(x,y2);
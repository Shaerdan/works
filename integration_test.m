%%%%% integration test %%%%%

clc;
clear all;
close all;

a=0.6;
f1 = @(x) log(1+a^2-2*a*cos(x));

f_intgral = integral(@(x) f1(x), -pi,pi);

theta = -pi:0.1:pi;
ff1 = f1(theta);
ff2 = 1+a^2-2*a*cos(theta);

plot(theta,ff1); hold on; plot(theta,ff2)
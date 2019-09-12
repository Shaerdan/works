clc;
clear all;
close all;

a= 1:1:2000;

b = a/pi - floor(a/pi);
plot(b)

y = sin(a);
y = y.^2;
%%%% contour of 2 d function:
clc;
clear all;
close all;

% f = @(x,y) x*sin(y);
x = linspace(0,1,100);
y = linspace(-pi,pi,100);
[X,Y] = meshgrid(x,y);

Z = (X.^2).*sin(Y);

contour(X,Y,Z,10,'ShowText','on')



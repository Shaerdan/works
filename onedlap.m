clc;
clear all;
close all;

global x1 x2 eps
x1 =0.9;
x2 =0.6;
eps=10^(-3.3);
h=0.0011;

d1 = 0.5*(x2^2-x1^2-x2+x1);
d2 = -0.5*((1-x2)^2 -(1-x1)^2);
d3 = 0.5*(x2^2-x1^2);


xx=0:h:1;

y1=0.5*((x2-xx).*sign(x2-xx)-(x1-xx).*sign(x1-xx));
figure(1)
plot(xx,y1)
test= mean(y1);
test2=mean(y1-d1)

y2 =  double((x2-xx).*heaviside(x2-xx)-(x1-xx).*heaviside(x1-xx));
figure(2)
plot(xx,y2)
test3=mean(y2);
test4=mean(y2-d3)
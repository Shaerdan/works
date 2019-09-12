clc;
clear all;

h=0.005;
a=0.5;
b=0.3;

x = 0.3:h:1;
y = 0.3:h:1;

[X,Y]=meshgrid(x,y);

f1 = -log(X.^2 + Y.^2);
df1 = del2(f1,h);

f2 = log((X-a).^2 + (Y-b).^2);
df2 = del2(f2,h);


k = 1/sqrt(a^2+b^2);
f3 = log((X-k*a).^2 + (Y-k*b).^2);
df3 = del2(f2,h);

ss=df2+df3;

surf(ss)
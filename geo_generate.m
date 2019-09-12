clc;
clear all;
close all;

r=1;               %radius of the disk
n1=100;             %boundary mesh number

h1=2*pi/n1;
theta=-pi:h1:pi;
[ x,y,xx,yy ] = mesh_circle( r,theta );
%check the plot
figure(1)
scatter(x,y);


%for non uniform mesh (reserved)

%Geometric parameters:
[ nx,ny,xm,ym,lx,ly,l ] = geom( x,y );
%check normal vectors:
figure(2)
quiver(x(1:lx),y(1:ly),nx,ny);



save('geo.mat')
function [ x,y,xx,yy ] = mesh_circle( r,theta )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

x = r*cos(theta);
y = r*sin(theta);
xx = x;
yy = y;

end


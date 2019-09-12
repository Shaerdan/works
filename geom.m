function [ nx,ny,xm,ym,lx,ly,l ] = geom( x,y )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
lx = length(x);
ly = length(y);
% for circles last point return to the first point
for i = 1:(lx-1)
    dy(i) = y(i+1)-y(i);
    dx(i) = x(i+1)-x(i);
end
    dy(ly) = dy(1);
    dx(lx) = dx(1);
    l = sqrt(dx.^2+dy.^2);
    nx = dy./l; 
    ny = -dx./l;
% Compute the midpoints:

for i = 1:lx-1
    xm(i) = (x(i+1)+x(i))/2;
    ym(i) = (y(i+1)+y(i))/2;
end

end


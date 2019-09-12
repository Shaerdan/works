function [ Se ] = fixloc_getSe(  n_e,n_source,xd,yd,x_plot,y_plot )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
       for i = 1:n_e
           for j = 1:n_source
               a = (xd(i)-x_plot(j))^2;
               b = (yd(i)-y_plot(j))^2;
               e = log(sqrt(a + b));
               Se(i,j) = (e/(2*pi));
           end
       end

end


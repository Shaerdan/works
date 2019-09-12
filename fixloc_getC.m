function [ C ] = fixloc_getC( n_mesh,n_source,xm,ym,x_plot,y_plot)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
for i = 1:n_mesh
    for j = 1 : n_source
        a = (xm(i) - x_plot(j))^2;
        b = (ym(i) - y_plot(j))^2;
        C(i,j) = -1/((2*pi)*sqrt(a + b));
    end
end

end


function [ output_args ] = fixloc_plot( S_plot1,S_plot2,x_plot,y_plot,...
                           x_fxloc,y_fxloc )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
% [X,Y,Z] = meshgrid(x_plot,y_plot,S);
% surf(X,Y,Z);
figure(1)
imagesc(x_plot,y_plot,S_plot1)

figure(2)
imagesc(x_plot,y_plot,S_plot2)

figure(3)
surf(x_fxloc,y_fxloc,S_plot1);

figure(4)
surf(x_fxloc,y_fxloc,S_plot2);

end


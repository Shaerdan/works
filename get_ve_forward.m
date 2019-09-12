function [ ve ] = get_ve_forward( xd,yd,S, x_plot, y_plot )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

for i = 1: length(S)
dist = sqrt((xd-x_plot(i)).^2 + (yd-y_plot(i)).^2);
ve(:,:,i) = (1/(2*pi))*S(i)*log(dist);
end
ve = sum(ve,3);

end


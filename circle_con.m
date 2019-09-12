function [ c,ceq ] = circle_con(x,n_source)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

S = x(1:n_source);
x_plot = x(n_source+1:2:end-1);
y_plot = x(n_source+2:2:end);
for i = 1:length(x_plot)
a(i) = (x_plot(i))^2 + (y_plot(i))^2 -1;
end

% for i = 1:length(x_plot)
% b(i) = 0.1 - ((x_plot(i))^2 + (y_plot(i))^2);
% end
c = a;
ceq = [];
end


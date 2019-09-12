function [ c,ceq ] = circle_con2(x,n_source)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


k=1;
for i = 1:n_source
    rs(i) = x(k+1);
a(i) = rs(i) -1;
k = k+3;
end

% for i = 1:length(x_plot)
% b(i) = 0.1 - ((x_plot(i))^2 + (y_plot(i))^2);
% end
c = a;
ceq = [];
end


function [phi] = phi_func(xx)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
global n_source S_real ths_real rs_real
for i = 1:n_source 
phi(:,i) = (S_real(i)/(2*pi))*log(sqrt(1+rs_real(i)^2 ...
    -2*rs_real(i)*cos(xx-ths_real(i))));
end

phi = 2*sum(phi,2)';
end


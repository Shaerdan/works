function [ phim, phi ] = boundary_src( r,theta,S,xi,yi )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Formulate BC
% Case 1: phi(x) known dphi(x) to be solved:
% Test 1: phi is constant
% phi = 1+zeros(1,lx-1);
% Test 2: phi is periodic, no source term: 
global alph
phi = r*cos(alph*theta);
% End Boundary Data
  temp = sqrt((r*cos(theta)-xi).^2 +...
      (r*sin(theta)-yi).^2);
  fsource = (S/(2*pi))*log(temp);
for i = 1:(length(theta)-1)  
phim(i) = (phi(i+1) + phi(i) - fsource(i+1) - fsource(i))/2;
end


end


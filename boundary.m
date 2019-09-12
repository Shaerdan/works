function [ phim, phi ] = boundary_phi( r,theta )
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
for i = 1:(length(theta)-1)
phim(i) = (phi(i+1) + phi(i))/2;
end


end


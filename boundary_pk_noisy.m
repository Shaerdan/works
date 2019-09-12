function [ pkm, pk ] = boundary_pk_noisy( nx,ny,x,y,xi,yi,S,snr )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Formulate BC
% Case 1: phi(x) known dphi(x) to be solved:
% Test 1: phi is constant
% phi = 1+zeros(1,lx-1);
% Test 2: phi is periodic, no source term: 

temp1 = (x-xi).^2 + (y-yi+eps(1)).^2;
temp2 = (x-xi).*nx + (y-yi).*ny;
gx = temp2./temp1;
pk = -((S+eps(2))/(2*pi))*gx ;
pk = awgn(pk,snr);
% End Boundary Data
for i = 1:(length(pk)-1)
pkm(i) = (pk(i+1) + pk(i))/2;
end


end


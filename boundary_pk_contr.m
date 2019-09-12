function [ pkm, pk ] = boundary_pk_contr( r,theta,nx,ny,x,y,x_source,S,sanity,dim )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Formulate BC
     global alph
if (sanity == 0)
    xi = x_source(1); 
    yi = x_source(2);
    temp1 = (x-xi).^2 + (y-yi).^2;
    temp2 = (x-xi).*nx + (y-yi).*ny;
    gx = temp2./temp1;
    pk= -(S/(2*pi))*gx;

end

if (sanity == 1)
pk = (alph*r^(alph-1))*cos(alph*theta);
end
    
% End Boundary Data
for i = 1:(length(pk)-1)
pkm(i) = (pk(i+1) + pk(i))/2;
end


end


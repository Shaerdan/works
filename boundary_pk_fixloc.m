function [ pkm, pk ] = boundary_pk_fixloc( nx,ny,x,y,x_plot,y_plot,S )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Formulate BC
    ll = length(x_plot);
    for i=1:ll
    temp1 = (x-x_plot(i)).^2 + (y-y_plot(i)).^2;
    temp2 = (x-x_plot(i)).*nx + (y-y_plot(i)).*ny;
    gx = temp2./temp1;
%     pk(k,:) = (1/(2*pi))-(S(k)/(2*pi))*gx;
    pk(i,:) = - (S(i)/(2*pi))*(1./(temp1).^(1/2));
    end
    pk =  sum(pk,1);
%         save('check.dat','gx');    
% End Boundary Data
for i = 1:(length(pk)-1)
pkm(i) = (pk(i+1) + pk(i))/2;
end


end




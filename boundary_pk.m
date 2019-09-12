function [ pkm, pk ] = boundary_pk( r,theta,n_source,nx,ny,x,y,S,x_plot,y_plot,sanity )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Formulate BC
     global alph
if (sanity == 0)
    for i=1:n_source
    temp1 = (x-x_plot(i)).^2 + (y-y_plot(i)).^2;
    temp2 = (1-x*x_plot(i)-y*y_plot(i));
%     pk(k,:) = (1/(2*pi))-(S(k)/(2*pi))*gx;
    pk(i,:) = - (S(i)/(2*pi))*(temp2./temp1);
    end
    pk =  sum(pk,1);
%         save('check.dat','gx');
end

if (sanity == 1)
% pk = (alph*r^(alph-1))*cos((alph+1)*theta);
pk = 2*(x.^2 - y.^2);
end

if (sanity == 2)
pk = (alph*1^(alph-1))*cos((alph)*theta);
end
    
% End Boundary Data
for i = 1:(length(pk)-1)
pkm(i) = (pk(i+1) + pk(i))/2;
end


end


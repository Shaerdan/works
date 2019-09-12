function [ Me,Ne ] = fixloc_getMeNe( n_e,xd,yd,nx,ny,l,xx,yy,lx,ly )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    for m=1:n_e
          F = integrals(l,nx,ny,xx,yy,xd(m),yd(m),lx,ly);
          Me(m,:) = F(2,:);
          Ne(m,:) = F(1,:);
    end

end


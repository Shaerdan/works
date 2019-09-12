function [ phim ] = Sol_BEM_Neumann( S,xi,yi,x,y,xx,yy,nx,ny,xm,ym,lx,...
                                 ly,l,r,theta,sanity,xc,yc,s1,s2,step,...
                                 )
                                    
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

   global flag1

    [ pkm, pk ] = boundary_pk( r,theta,nx,ny,x,y,xi,yi,S,sanity );
    
      for i = 1: lx-1     % loop over mid point
      [F1,delta] = integrals(l,nx,ny,xx,yy,xm(i),ym(i),lx,ly);
        for j = 1:lx-1   % loop over boundary mesh
            delt = (i==j)*1;
            A(i,j) = (F1(2,j) - (1/2)*delt);
            B(i,j) = pkm(j)*F1(1,j);
        end
      end
  b = sum(B,2);
%   restart = 10;
%   tol = 10^(-8);
%   maxit = 100;
%   z = gmres(A,b,restart,tol,maxit);
shift = 0;
tol = 10^(-12);
[z,flag1] = cgls(A,b,shift,tol);
  phim = z';
  
  phim_measure = phim(1:step:end);
  
  

%     for m=1:s1-1
%         for n = 1:s2
%                 F = integrals(l,nx,ny,xx,yy,xc(m,n),yc(m,n),lx,ly);
%                 soln(m,n) =  F(2,:)*phim' - F(1,:)*pkm';
%         end
%     end
  
  
end


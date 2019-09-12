function [ F1,A,P ] = DtN( xx,yy,...
                              nx,ny,xm,ym,lx,...
                              ly,l,...
                              gmres_restart,gmres_tol,gmres_maxit)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

      for i = 1: lx-1     % loop over mid point
      [F1,delta] = integrals(l,nx,ny,xx,yy,xm(i),ym(i),lx,ly);
        for j = 1:lx-1   % loop over boundary mesh
            delt = (i==j)*1;
            A(i,j) = -F1(1,j);
            B(i,j) = (-F1(2,j)+0.5*delt);
            P(i,j) = (-F1(2,j)+0.5*delt);
        end
      end
%   b = sum(B,2);
  cond(P)
  cond(A)
  
%   z = gmres(L,g',gmres_restart,gmres_tol,gmres_maxit);
%   ff = z';
% shift = 0;
% tol = 10^(-10);
% % cond(A)
% [z,flag1] = cgls(A,b,shift,tol);
%   gg = z';    
end


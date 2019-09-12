function [ ff,F1,A,P ] = NtD_miu2( g,xx,yy,...
                                   nx,ny,xm,ym,lx,...
                                   ly,l,reg_miu,K1,...
                                   gmres_restart,gmres_tol,gmres_maxit)
                                    
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

      for i = 1: lx-1     % loop over mid point
      [F1,delta] = integrals(l,nx,ny,xx,yy,xm(i),ym(i),lx,ly);
        for j = 1:lx-1   % loop over boundary mesh
            delt = (i==j)*1;
            A(i,j) = (F1(2,j) - (1/2)*delt);
            B(i,j) = g(j)*F1(1,j);
            P(i,j) = F1(1,j);
        end
      end
%       b = sum(B,2);       
       L = inv(P)*A;
       b = L'*g'+ reg_miu*L'*K1;
       S = L'*L;

  z = gmres(S,b,gmres_restart,gmres_tol,gmres_maxit);
  ff = z';
% shift = 0;
% tol = 10^(-16);
% % cond(A)
% [z,flag1] = cgls(S,b,shift,tol);
%   ff = z';
%   flag1
end


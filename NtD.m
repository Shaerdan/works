function [ ff,F1,A,P ] = NtD( g,xx,yy,...
                              nx,ny,xm,ym,lx,...
                              ly,l,alph_tik)
                                    
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
      b = sum(B,2);
%       b = P*g';


%   restart = 25;
%   tol = 10^(-16);
%   maxit = 500;
%   z = gmres(A,b,restart,tol,maxit);

% shift = 0;
% tol = 10^(-16);
% % cond(A)
% [z,flag1] = cgls(A,b,shift,tol);
%   ff = z';

%   flag1
%   if(flag1 ~= 1)
%       flag1
%   end
%%%%%Tichnov Regularisation:
[Ua S_a Va] = svd(A);
Sa_diag = diag(S_a);
Dap = diag(Sa_diag./(Sa_diag.^2 + alph_tik^2));
Ar_inv = Va*Dap*Ua';
z = Ar_inv*b;
ff = z';  

% nm = length(b);
% %%%%% Deflation with Tichnov %%%%
% GAMMA_reg = A'*A + alph_tik*eye(nm,nm) +...
%      alph_deflate*(eye(nm,1)*eye(1,nm));
% ff =  GAMMA_reg\(A'*b);
% ff = ff';
  
end


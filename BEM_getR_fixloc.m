function [ R ] = BEM_getR_fixloc( S,x_plot,y_plot,x,y,xx,yy,nx,...
                         ny,xm,ym,lx,ly,l,sanity,xd,yd,...
                         f_measure_noisy,lambda_inv,alph_tik)
                                    
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%     S = source_sol(1:n_source);
%     x_source = source_sol(n_source+1:end);
    
    [ g, pk ] = boundary_pk_fixloc( nx,ny,x,y,x_plot,y_plot,S);
    
      for i = 1: lx-1     % loop over mid point
      [F1,delta] = integrals(l,nx,ny,xx,yy,xm(i),ym(i),lx,ly);
        for j = 1:lx-1   % loop over boundary mesh
            delt = (i==j)*1;
            A(i,j) = (F1(2,j) - (1/2)*delt);
            P(i,j) = F1(1,j);
            B(i,j) = g(j)*F1(1,j);
        end
      end
      b = sum(B,2);
%       b = P*g';

%   restart = 10;
%   tol = 10^(-8);
%   maxit = 100;
%   z = gmres(A,b,restart,tol,maxit);

% shift = 0;
% tol = 10^(-16);
% [z,flag1] = cgls(A,b,shift,tol);
% ff = z';  

% if (flag1 ~= 1)
%     flag1
% end

%%%%%Tichnov Regularisation:
[Ua, S_a, Va] = svd(A);
Sa_diag = diag(S_a);
Dap = diag(Sa_diag./(Sa_diag.^2 + alph_tik));
Ar_inv = Va*Dap*Ua';
z = Ar_inv*b;
ff = z';

  %%% Get the BEM solution at the measured data piont location
%   f_data = zeros(length(xd),1);
    n_e = length(xd);
    n_source = length(S);
    for m=1:n_e
          F = integrals(l,nx,ny,xx,yy,xd(m),yd(m),lx,ly);
          f_data(m) =  2*(F(2,:)*ff' - F(1,:)*g');
    end
    [ Se ] = fixloc_getSe(  n_e,n_source,xd,yd,x_plot,y_plot );
    ve = Se*S;
    f_data = f_data + ve';
%       f_data = 2*(M*ff'-N*g')' ;

%     R = norm(phim_data-phim_measure_noisy,2);
    % Tichnov Regularisation:
    R = norm(f_data - f_measure_noisy,2)^2 ...
        + lambda_inv*norm(f_data,2)^2;
    
  
end


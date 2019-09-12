function [ R ] = BEM_getR_Neumann_miu( source_sol,n_source,n_loc,x,y,xx,yy,nx,ny,xm,ym,lx,...
                                ly,l,r,theta,sanity,xd,yd,...
                                f_measure_noisy,dim,lambda_inv,reg_miu,K1,...
                                gmres_restart,gmres_tol,gmres_maxit)
                                    
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    S = source_sol(1:n_source);
    x_source = source_sol(n_source+1:end);
    
    [ g, pk ] = boundary_pk( r,theta,n_source,nx,ny,x,y,S,x_plot,y_plot,sanity,dim );
    
%%%%%Miu regularisation:
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
  %%% Get the BEM solution at the measured data piont location
    for m=1:length(xd)
          F = integrals(l,nx,ny,xx,yy,xd(m),yd(m),lx,ly);
          f_data(m) =  2*(F(2,:)*ff' - F(1,:)*g');
    end
%       f_data = 2*(M*ff'-N*g')' ;

%     R = norm(phim_data-phim_measure_noisy,2);
    % Tichnov Regularisation:
    R = (norm(f_data-f_measure_noisy,2))^2 ...
        + lambda_inv*(norm(f_data,2));
    
  
end


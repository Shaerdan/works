function [ R ] = BEM_getR_Dirichlet( phim,source_sol,n_source,n_loc,...
                                x,y,xx,yy,nx,ny,xm,ym,lx,...
                                ly,l,r,theta,sanity,xd,yd,...
                                phim_measure_noisy,dim,lambda1,lambda2)
                                    
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    S = source_sol(1:n_source);
    x_source = source_sol(n_source+1:end);
[ pkm,F1 ] = Sol_Poisson_BEM_Dirichlet(phim,...
                      x_source,S,x,y,xx,yy,nx,ny,xm,ym,lx,...
                      ly,l,r,theta,sanity,dim);
  %%% Get the BEM solution at the measured data piont location
    for m=1:length(xd)
                F = integrals(l,nx,ny,xx,yy,xd(m),yd(m),lx,ly);
                phim_data(m) =  2*(F(2,:)*phim' - F(1,:)*pkm');
    end
  
%     R = norm(phim_data-phim_measure_noisy,2);
    % Tichnov Regularisation:
    R = lambda1*norm(phim_data-phim_measure_noisy,2)...
        + lambda2*norm(phim_data,2);
    
  
end


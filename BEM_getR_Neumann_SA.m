function [ R ] = BEM_getR_Neumann_SA( source_sol,n_source,n_loc,x,y,xx,yy,nx,ny,xm,ym,lx,...
                                ly,l,r,theta,sanity,xd,yd,...
                                phim_measure_noisy,dim,lambda1,lambda2,weight1,...
                                weight2)
                                    
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    S = source_sol(1:n_source);
    x_source = source_sol(n_source+1:end);

    [ pkm, pk ] = boundary_pk( r,theta,nx,ny,x,y,x_source,S,sanity,dim );
    
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
tol = 10^(-10);
[z,flag1] = cgls(A,b,shift,tol);
phim = z';  
if (flag1 ~= 1)
    flag1
end
  %%% Get the BEM solution at the measured data piont location
    for m=1:length(xd)
                F = integrals(l,nx,ny,xx,yy,xd(m),yd(m),lx,ly);
                phim_data(m) =  2*(F(2,:)*phim' - F(1,:)*pkm');
    end
  
%     R = norm(phim_data-phim_measure_noisy,2);
    % Tichnov Regularisation with constraint:
    Aa = [(ones(n_source,1))' (zeros(2*n_source,1))'];    % Make it flexible for more than 2 sources
    ba = [0];
    f_c = @(x) circle_con(x);    
    penalty1 = Aa*source_sol - ba;
    penalty2 = max(f_c(source_sol),0)
    R = lambda1*norm(phim_data-phim_measure_noisy,2)...
        + lambda2*norm(phim_data,2)+weight1*norm(penalty1,2)...
        +weight2*sum(penalty2);
    
  
end




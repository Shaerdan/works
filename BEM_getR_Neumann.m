function [ R ] = BEM_getR_Neumann( source_sol,n_source,n_loc,...
                                x,y,xx,yy,nx,...
                                ny,xm,ym,lx,ly,l,r,theta,sanity,xd,yd,...
                                f_measure_noisy,dim,lambda_inv,alph_tik);
                                    
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

S = source_sol(1:n_source);
x_source = source_sol(n_source+1:2:end-1);
y_source = source_sol(n_source+2:2:end);

    [ g, pk ] = boundary_pk( r,theta,n_source,nx,ny,x,y,S,x_source,y_source,sanity );
    
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

%%%%%Tichnov Regularisation:
[Ua S_a Va] = svd(A);
Sa_diag = diag(S_a);
Dap = diag(Sa_diag./(Sa_diag.^2 + alph_tik));
Ar_inv = Va*Dap*Ua';
z = Ar_inv*b;
ff = z';

  %%% Get the BEM solution at the measured data piont location
    for m=1:length(xd)
          F = integrals(l,nx,ny,xx,yy,xd(m),yd(m),lx,ly);
          f_data(m) =  2*(F(2,:)*ff' - F(1,:)*g');
    end
        [ ve ] = get_ve( xd,yd,S, x_source, y_source );
  
        f_data = f_data + ve;
%       f_data = 2*(M*ff'-N*g')' ;

%     R = norm(phim_data-phim_measure_noisy,2);
    % Tichnov Regularisation:
    R = norm(f_data-f_measure_noisy,2)^2 ...
        + lambda_inv*norm(f_data,2)^2;
    
    % Mumford-Shah:
%     for i = 1:n_source
%         rs(i) = sqrt(x_source(i)^2 + y_source(i)^2);
% fun = @(xx) (sin(xx-ths(i))./(1+rs(i)^2-2*rs(i)*cos(xx - ths(i)))).^2;
% reglr(i) = (1/2*pi)*S(i)*rs(i)*integral(fun,-pi,pi);
% end
% reglr = sum(reglr)
% 
%     R = norm(f_data-f_measure_noisy,2)^2 ...
%         + lambda_inv*reglr;
    
    
  
end


clc;
clear all;
close all;

% Control parameters:
sanity = 0;        %% Dummy, change it
type = 'N';     % 1.Boundary type switch,  D--> Dirichlet, N --> Neumann.  
constrt = 1;
method = 1;    % 1 for fmincon, 2 for simulated annealing, 3 for hybrid
lambda_inv = 0.1; % Tichnov Regularisation for the inverse problem
%%%% Use geo_generate routine to generate artificial data
%%%% Load geometry 
load('geo.mat')
load('geo2.mat')

%%%% Use data_generate routine to generate artificial data
%%%% Load data
load('datafixloc.mat')

                                    
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
    for m=1:length(xd)
          F = integrals(l,nx,ny,xx,yy,xd(m),yd(m),lx,ly);
          Me(m,:) = F(2,:);
          Ne(m,:) = F(1,:);
          f_data(m) =  2*(F(2,:)*ff' - F(1,:)*g');
    end
    f_data2 = 2*(Me*ff'-Ne*g')' ;

%     R = norm(phim_data-phim_measure_noisy,2);
    % Tichnov Regularisation:
    R = norm(f_data - f_measure_noisy,2)^2 ...
        + lambda_inv*norm(f_data,2)^2;
    


% size1 = length(S)/2;
% size2 = length(S) - size1;
% img1 = reshape(S,size1,size2);
% img2 = reshape(x_invsoln,size1,size2);
% figure(6)
% trisurf(tri,x_plot,y_plot,S);
% caxis([-1 1])
% figure(7)
% trisurf(tri,x_plot,y_plot,x_invsoln);
% caxis([-1 1])


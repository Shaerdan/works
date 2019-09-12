clc;
clear all;
close all;

% Control parameters:
global alph 
alph = 2;          %dimention parameter for periodic bc in sanity test
levels = 40;       %error contour levels
% Soure term parameters, take S = 0 for no source term case:
load('data.mat');
% Switches:
% 1.Boundary type switch, pk_phi = 0 --> Dirichlet, pk_phi = 1 --> Neumann:    
pk_phi = 1;

% 2.Sanity Test Switch, if sanity == 1, perform u = u(r) = r^(alpha)*cos(alpha*theta)
% test, Dirichlet boundary condition is u = u(r_bc), in Neumann bc is
% du(r)/dr at r_bc;
sanity = 0;

% 3. Error Contour Swith: %%%%% MULTI SOURCE CONTOUR NOT WORKING YET
contr = 1;   

% Get geometrical data:
load('geo.mat')


if ( pk_phi ==0 )     % Dirichlet BC Problem Case
[phim,phi] = boundary_phi(r,theta);

  for i = 1: lx-1     % loop over mid point
    [F1,delta] = integrals(l,nx,ny,xx,yy,xm(i),ym(i),lx,ly);
      for j = 1:lx-1   % loop over mesh
          delt = (i==j)*1;
          A(i,j) = -F1(1,j);
          B(i,j) = phim(j)*(-F1(2,j) + (1/2)*delt);
      end
  end
  b = sum(B,2);

  z = gmres(A,b);
  pkm = z';
  
end

if ( pk_phi ==1 )      % Neumann Boundary Condition Case
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
tol = 10^(-12);
[z,flag1] = cgls(A,b,shift,tol);
  phim = z';
end


% z = linsolve(A,b);
% Other linear solvers:


% Solve the internel solution loop over all internel verteces;

% % For radiate grids (from center):
hh = r/10;
rho = 0:hh:r;
[th, rr] = meshgrid(theta, rho(1:end-1));
xc = rr.*cos(th);
yc = rr.*sin(th);


s1 = length(rho);
s2 = length(theta);

% % For triangular grids (NOT FINISHED):
% fd=@(p) sqrt(sum(p.^2,2))-1;
% [p,t]=distmesh2d(fd,@huniform,0.05,[-0.5,-0.5;0.5,0.5],[]);
% xc = p(:,1);
% yc = p(:,2);


    for m=1:s1-1
        for n = 1:s2
                F = integrals(l,nx,ny,xx,yy,xc(m,n),yc(m,n),lx,ly);
                soln(m,n) =  (F(2,:)*phim' - F(1,:)*pkm');
        end
    end
% Extract boundary data from boundary average data:
% M = zeros(s2-1,s2);
% for i = 1:(s2-1)
%      M(i,i)   = 1/2;
%      M(i,i+1) = 1/2;
% end
% shift = 0;
% tol = 10^(-16);
% maxit = 500;
% [phi_bc,flag2] = cgls(M,phim',shift,tol,maxit);
%     for ss=1:length(x)
%           F = integrals(l,nx,ny,xx,yy,x(ss),y(ss),lx,ly);
%           phi_bc(ss) =  F(2,:)*phim' - F(1,:)*pkm';
%     end
% soln(end+1,:) = phi_bc';
    ll = length(x_source);
    k=1;
    if (sanity == 0)
    for i = 1:dim:ll
    xi = x_source(i); yi = x_source(i+1);
    N(:,:,k) = (S(k)/(2*pi))*log(sqrt((xc-xi).^2 +(yc-yi).^2)); 
    k=k+1;
    end
    N = sum(N,3);
    u = N + soln; 
    end

    if (sanity == 1)
        u = soln;
       % Analytical Soln:
    soln_a = (rr.^alph).*cos(alph*th);
        figure(4)
       surf(xc,yc,soln_a)
       zlim([-1 1])
       
       % Error matrix:
       err = abs(soln_a - soln);       
       err_rel = err./soln_a;
    end
    
    figure(5)
       surf(xc,yc,u)
%        zlim([-1 1])


       condA = cond(A);




    
    
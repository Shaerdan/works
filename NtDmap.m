clc;
clear all;
close all;

% Control parameters:
global alph 
alph = 2;          %dimention parameter for periodic bc in sanity test
levels = 40;       %error contour levels
% Get data:
load('data.mat');
load('geo.mat')
% pk_phi = 1;
lambda_NtD = 0;
% 2.Sanity Test Switch, if sanity == 1, perform u = u(r) = r^(alpha)*cos(alpha*theta)
% test, Dirichlet boundary condition is u = u(r_bc), in Neumann bc is
% du(r)/dr at r_bc;
sanity = 0;

% 3. Error Contour Swith: %%%%% MULTI SOURCE CONTOUR NOT WORKING YET
contr = 1;   

[ pkm, pk ] = boundary_pk( r,theta,nx,ny,x,y,x_source,S,sanity,dim );
phim_0 = zeros(size(pkm)); 
% phim_0 = g;
pkm_match = pkm(1:stp:end);

 f = @(phim)R_NtD( pkm,pkm_match,phim,phim_measure_noisy,...
                      x_source,S,x,y,xx,yy,nx,ny,xm,ym,lx,...
                      ly,l,r,theta,sanity,dim,stp,lambda_NtD );   
% options = optimset('PlotFcns',@optimplotfval);   
% [g_0,fval,exitflag,output] = fminsearch(f,phim_0,options);

options = optimoptions(@fminunc,'PlotFcns',@optimplotfval,'Algorithm','quasi-newton');
[g,fval,exitflag,output] = fminunc(f,phim_0,options);
% 
% options = optimset('PlotFcns',@optimplotfval,'MaxIter',5000);   
% [g,fval,exitflag,output] = fminsearch(f,g_0,options);
% figure(4)
% plot(phim_0);
figure(1)
plot(phim_data); hold on; plot(g);
errrel = abs(g-phim_data)./abs(phim_data);
save('NtD','g');

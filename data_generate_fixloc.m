clc;
close all;
clear all;

global alph snr
snr = 'no';          %signal to noise ratio for artificial noisy boundary%
alph = 2;          %dimention parameter for periodic bc in sanity test
alph_tik = 0.001;   % Tichnov Regularisation for NtD map
data_pts = 15;
data_perccet=1/data_pts;      
sanity = 0;
lump = 10;
full_partial = 0;   % 1 for full data, 0 for partial;
manuel = 2;         % 0 for full random, 1 for full manuel, 2 for symmetric
                    % 3 for unsymmetric
reg = 0;            % 0 for Tichnov NtD, 1 for miu regularised Ntd
                    % 2 for squared miu with kernel
reg_miu =-0.1;
random = 2;

%%control pamameters for gmres solver:
gmres_restart = 25;
gmres_tol = 10^(-10); 
gmres_maxit = 1000;  
load('geo.mat')
% load('geo2.mat')
%     n_source = length(x_plot);

    if (random ~= 3)
    load('geo2.mat');
    n_source = length(x_plot);
    end
    if (random == 1)
    postv = 0.5;
    
    n_sp = (postv*n_source); 
    n_sm = n_source - n_sp;
    S1 = randn(n_sp,1);
    S1 = S1/sum(S1);
%     S2 = -rand(n_sm,1);
    S2 = -S1;
    S = [S1;S2];

    elseif (random == 2)
    S = zeros(n_source,1);
    S(1:lump) = 1;
    S(end-lump:end) = -1;
    elseif (random == 3)
       n_sp = 3;   
       n_sm = 3; 
        S = [1 1 -1 1 -1 -1]';
       
       n_source = length(S);
       ang_stp = (2*pi)/n_source;
       th = 0:ang_stp:2*pi;
       rr = 0.85;

        for i = 1:n_source
         x_plot(i) = rr*cos(th(i));
         y_plot(i) = rr*sin(th(i));
        end 
    end   
    
    

[ g, pk ] = boundary_pk_fixloc(nx,ny,x,y,x_plot,y_plot,S);
C = fixloc_getC( n1,n_source,x,y,x_plot,y_plot );
g_C = (1/(2*pi))*C*S;

   if(reg == 0)
    [ ff,F1,A,P ] = NtD( g,xx,yy,...
                         nx,ny,xm,ym,lx,...
                         ly,l,alph_tik);

    elseif (reg == 1)
        Ka = tril(ones(n,n),-1);
        K1 = ones(n,1);
        [ ff,F1,A,P ] = NtD_miu2( g,xx,yy,...
                                   nx,ny,xm,ym,lx,...
                                   ly,l,reg_miu,K1,...
                                   gmres_restart,gmres_tol,gmres_maxit);
    elseif (reg == 2)
        Ka = tril(ones(n,n),-1);
        K1 = ones(n,1);
        [ ff,F1,A,P ] = NtD_miu( g,xx,yy,...
                                   nx,ny,xm,ym,lx,...
                                   ly,l,reg_miu,Ka,...
                                   gmres_restart,gmres_tol,gmres_maxit);
                                    
    end
    stp = data_perccet*n1;
    xd = xm(1:stp:end);
    yd = ym(1:stp:end);
    n_e = length(xd);
    [ Se ] = fixloc_getSe(  n_e,n_source,xd,yd,x_plot,y_plot );
    
    ve = Se*S;
%     for m=1:length(xm)                               
    for m=1:length(xd)
          F = integrals(l,nx,ny,xx,yy,xd(m),yd(m),lx,ly);
          f_data(m) =  2.0*(F(2,:)*ff' - F(1,:)*g');
          M(m,:) = F(2,:);
          N(m,:) = F(1,:);
    end    
        
    f_data = f_data + ve';
    
    if (full_partial == 1)
    for m=1:length(xm)                               
%     for m=1:length(xd)
          F = integrals(l,nx,ny,xx,yy,xm(m),ym(m),lx,ly);
          f_data(m) =  2*(F(2,:)*g' - F(1,:)*ff');
    end
    end
   
    if (strcmp(snr,'no'))
    f_measure_noisy = f_data;
    else   
    f_measure_noisy = awgn(f_data,snr);
    norm(f_measure_noisy,2)
    end
    
    figure(2)
    plot(f_data); hold on; plot(f_measure_noisy);

    save('datafixloc.mat')
                                       
                                       
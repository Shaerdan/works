%%%% Create artificial measurement:
%%%% Run forward problem once to get boundary data from a set 
clc;
close all;
clear all;

global alph snr flag1 
snr = 'no';          %signal to noise ratio for artificial noisy boundary%
alph = 2;          %dimention parameter for periodic bc in sanity test
alph_tik = 0.0001;   % Tichnov Regularisation for NtD map
alph_deflate = 10^8; % deflation of NtD map

lambda_inv = 0.001; % Tichnov Regularisation for the inverse problem
data_pts = 100;
sanity = 0;
data_perccet=1/data_pts;      
full_partial = 0;   % 1 for full data, 0 for partial;
manuel = 1;         % 0 for full random, 1 for full manuel, 2 for symmetric
                    % 3 for unsymmetric
if (manuel == 1|2)
n_sp = 3;           % number of points: n_sp*2         
rr = 0.7;
end
if (manuel == 3)
rr = 0.3;         % radius range with normal distribution
sigma = 0.0;        % deviation of the radius range
n_sp = 2;          % positive charges
n_sm = 2;          % negative charges
fan_range = 2*pi;  % scattered range angle 
end
reg = 0;            % 0 for Tichnov NtD, 1 for miu regularised Ntd
                    % 2 for squared miu with kernel
reg_miu =0;

%%control pamameters for gmres solver:
gmres_restart = 25;
gmres_tol = 10^(-10); 
gmres_maxit = 1000;
%%%% Use geo_generate routine to generate artificial data
%%%% Load geometry 


load('geo.mat');
    dim = 2;
   
   if(manuel == 1)
    S = [1 -1];
    n_source = length(S);
    ang1 = -pi; small_ang = 0.05;
    ang2 = -pi +small_ang;
    fan_range = ang2 - ang1;
    ang_stp = (fan_range)/(n_source-1);
    th = ang1:ang_stp:ang2;
%     th = [th 2*pi];
    rr = 0.7;
    
        for i = 1:n_source
         x_plot(i) = rr*cos(th(i));
         y_plot(i) = rr*sin(th(i));
        end
        
        for i = 1:n_source-1
        dis(i) = (x_plot(i+1)-x_plot(i))^2 +(y_plot(i+1)-y_plot(i))^2;
        end
        dis = sqrt(sum(dis));
%     plot(x,y); hold on;
   end 
     k     =   1;
     n_loc        = 2*length(x_plot);

 %%%% Forward BEM method to solve boundary data phim 
    [ g, gk ] = boundary_pk( r,theta,n_source,nx,ny,x,y,S,x_plot,y_plot,sanity );

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
    
    [ ve ] = get_ve( xd,yd,S, x_plot, y_plot );
%     for m=1:length(xm)                               
    for m=1:length(xd)
          F = integrals(l,nx,ny,xx,yy,xd(m),yd(m),lx,ly);
          f_data(m) =  2*(F(2,:)*ff' - F(1,:)*g');
          M(m,:) = F(2,:);
          N(m,:) = F(1,:);
    end  
    
        figure(4)
    plot(f_data);hold on; plot(ve)
    
    f_data = f_data + ve;
    
   
    if (strcmp(snr,'no'))
    f_measure_noisy = f_data;
    else   
    f_measure_noisy = awgn(f_data,snr);
    norm(f_measure_noisy,2)
    end
    
    figure(1)
    subplot(2,2,1)
    plot(x,y); hold on;
    scatter(xd,yd,'*')
    k = 1;
    for i = 1:n_source
        scatter(x_plot(i),y_plot(i));
        text1 = S(k);
        text1 = num2str(text1);
        text1 = cellstr(text1);
        k = k + 1;
        text_fin =text(x_plot(i)+0.05,y_plot(i)+0.05,text1);
        fontsize = text_fin.FontSize;
        text_fin.FontSize = 4;
        hold on;
    end
    subplot(2,2,2)
    plot(f_data,'*'); hold on; plot(f_measure_noisy);
    
    figure(2)
    theta_sample=-pi:(2*pi/(m-1)):pi;
    plot(theta_sample,f_data,'-*');
    
    figure(3)
        plot(x,y); hold on;
    k = 1;
    for i = 1:n_source
        scatter(x_plot(i),y_plot(i));
        text1 = S(k);
        text1 = num2str(text1);
        text1 = cellstr(text1);
        k = k + 1;
        text_fin =text(x_plot(i)+0.05,y_plot(i)+0.05,text1);
        fontsize = text_fin.FontSize;
        text_fin.FontSize = 4;
        hold on;
    end


    
    if (reg == 0)
    save('data.mat','alph_tik','f_measure_noisy','f_data','x_plot','y_plot','S',...
        'n_source','dim','xd','yd','M','N','A','P','reg','gmres_restart',...
        'gmres_tol','gmres_maxit','lambda_inv','n_loc');    
    else
    save('data.mat','alph_tik','f_measure_noisy','f_data','x_plot','y_plot','S',...
        'n_source','dim','xd','yd','M','N','A','P','reg','reg_miu',...
        'Ka','K1','gmres_restart',...
        'gmres_tol','gmres_maxit','lambda_inv','n_loc')
    end
    % Get Data Variance Covariance matrix
%     Sig = (1/length(f_measure_noisy))*f_measure_noisy'*f_measure_noisy;
%     EigSig = eig(Sig);
    %%%% End of solving boundary data using BEM  
    
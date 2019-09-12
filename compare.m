    clc;
    clear all;
    close all;
    
    load('data.mat');
    load('geo.mat')
    load('inv_soln.mat');

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
     alph = 2;          %dimention parameter for periodic bc in sanity test
     alph_tik = 0.0001;   % Tichnov Regularisation for NtD map
     lambda_inv = 0.0; % Tichnov Regularisation for the inverse problem
     data_perccet=0.1;      
     sanity = 0;
    full_partial = 0;   % 1 for full data, 0 for partial;
    manuel = 2;      
%     n_locop        = length(x_source);

 %%%% Forward BEM method to solve boundary data phim 
    [ g, pk ] = boundary_pk( r,theta,nx,ny,x,y,loc_inv,S_inv,sanity,dim );

    if(reg == 0)
    [ ff,F1,A,P ] = NtD( g,loc_inv,S_inv,x,y,xx,yy,...
                                   nx,ny,xm,ym,lx,...
                                   ly,l,r,theta,sanity,dim,alph_tik);

    elseif (reg == 1)
        Ka = tril(ones(n,n),-1);
        K1 = ones(n,1);
        [ ff,F1,A,P ] = NtD_miu2( g,loc_inv,S_inv,x,y,xx,yy,...
                                   nx,ny,xm,ym,lx,...
                                   ly,l,r,theta,sanity,dim,reg_miu,Ka,K1,...
                                   gmres_restart,gmres_tol,gmres_maxit);
    elseif (reg == 2)
        Ka = tril(ones(n,n),-1);
        K1 = ones(n,1);
        [ ff,F1,A,P ] = NtD_miu( g,loc_inv,S_inv,x,y,xx,yy,...
                                   nx,ny,xm,ym,lx,...
                                   ly,l,r,theta,sanity,dim,reg_miu,Ka,K1,...
                                   gmres_restart,gmres_tol,gmres_maxit);
                                    
    end
    stp = data_perccet*n;
    xd = xm(1:stp:end);
    yd = ym(1:stp:end);
%     for m=1:length(xm)                               
    for m=1:length(xd)
          F = integrals(l,nx,ny,xx,yy,xd(m),yd(m),lx,ly);
          f_data_reconstr(m) =  2.0*(F(2,:)*ff' - F(1,:)*g');
          M(m,:) = F(2,:);
          N(m,:) = F(1,:);
    end     
    if (full_partial == 1)
    for m=1:length(xm)                               
%     for m=1:length(xd)
          F = integrals(l,nx,ny,xx,yy,xm(m),ym(m),lx,ly);
          f_data_reconstr(m) =  2*(F(2,:)*g' - F(1,:)*ff');
    end
    end
    err_data = f_data - f_data_reconstr;
    rel_err_data = err_data./(abs(f_data));
    figure(1)
    plot(f_data); hold on; plot(f_data_reconstr);
    figure(2)
    plot(err_data);
    figure(3)
    plot(rel_err_data);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
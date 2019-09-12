%%%% Create artificial measurement:
%%%% Run forward problem once to get boundary data from a set 
clc;
close all;
clear all;

global alph snr flag1 
snr = 50000000;          %signal to noise ratio for artificial noisy boundary%
alph = 2;          %dimention parameter for periodic bc in sanity test
alph_tik = 0.0001;   % Tichnov Regularisation for NtD map
lambda_inv = 0.0; % Tichnov Regularisation for the inverse problem
data_perccet=0.1;      
sanity = 0;
full_partial = 0;   % 1 for full data, 0 for partial;

manuel = 2;         % 0 for full random, 1 for full manuel, 2 for symmetric
                    % 3 for unsymmetric


reg = 0;            % 0 for Tichnov NtD, 1 for miu regularised Ntd
                    % 2 for squared miu with kernel
reg_miu =-0.1;

%%control pamameters for gmres solver:
gmres_restart = 25;
gmres_tol = 10^(-10); 
gmres_maxit = 1000;
%%%% Use geo_generate routine to generate artificial data
%%%% Load geometry 
load('geo.mat');
    dim = 2;
   

       rstr = 0.1;
       hh = 0.4;
       rend = 0.9;
    
       
       n_sp = 2;           % number of points: n_sp*2  
    
       n_sm = 2;          % negative charges
  
       S1 = rand(n_sp,1); S2 = -S1;
       S= [S1;S2];
       n_source = length(S);
       ang_stp = (2*pi)/n_source;
       th1 = 0:ang_stp:2*pi;
       
        stp = data_perccet*n;
        xd = xm(1:stp:end);
        yd = ym(1:stp:end);
       
        sigma = 0.0;        % deviation of the radius range
        fan_range = 0.5*pi;  % scattered range angle
        
        
       kk = 1;
 for rr = rstr:hh:rend
%         rr = 0.15;
       
        for i = 1:n_source
         x_source1(i) = rr*cos(th1(i));
         x_source2(i) = rr*sin(th1(i));
        end 
                k=1;
%     plot(x,y); hold on;
       for i = 1:dim:2*n_source-1
         x_source(i) = x_source1(k);
         x_source(i+1) = x_source2(k);
         k = k+1;
%     scatter( x_source(i), x_source(i+1)); hold on;
       end
     n_loc        = length(x_source);

 %%%% Forward BEM method to solve boundary data phim 
    [ g, pk ] = boundary_pk( r,theta,nx,ny,x,y,x_source,S,sanity,dim );

    if(reg == 0)
    [ ff,F1,A,P ] = NtD( g,x_source,S,x,y,xx,yy,...
                                   nx,ny,xm,ym,lx,...
                                   ly,l,r,theta,sanity,dim,alph_tik);

    elseif (reg == 1)
        Ka = tril(ones(n,n),-1);
        K1 = ones(n,1);
        [ ff,F1,A,P ] = NtD_miu2( g,x_source,S,x,y,xx,yy,...
                                   nx,ny,xm,ym,lx,...
                                   ly,l,r,theta,sanity,dim,reg_miu,Ka,K1,...
                                   gmres_restart,gmres_tol,gmres_maxit);
    elseif (reg == 2)
        Ka = tril(ones(n,n),-1);
        K1 = ones(n,1);
        [ ff,F1,A,P ] = NtD_miu( g,x_source,S,x,y,xx,yy,...
                                   nx,ny,xm,ym,lx,...
                                   ly,l,r,theta,sanity,dim,reg_miu,Ka,K1,...
                                   gmres_restart,gmres_tol,gmres_maxit);
                                    
    end

%     for m=1:length(xm)                               
    for m=1:length(xd)
          F = integrals(l,nx,ny,xx,yy,xd(m),yd(m),lx,ly);
          f_data(m) =  2.0*(F(2,:)*ff' - F(1,:)*g');
          M(m,:) = F(2,:);
          N(m,:) = F(1,:);
    end     
    if (full_partial == 1)
    for m=1:length(xm)                               
%     for m=1:length(xd)
          F = integrals(l,nx,ny,xx,yy,xm(m),ym(m),lx,ly);
          f_data(m) =  2*(F(2,:)*g' - F(1,:)*ff');
    end
    end
    
    f_measure_noisy = awgn(f_data,snr);
    norm(f_measure_noisy,2)
    
    figure(1)
    plot(x,y); hold on;
    for i = 1:dim:2*n_source
        scatter(x_source(i),x_source(i+1));
        hold on;
    end
    figure(2)
    plot(f_data); hold on; plot(f_measure_noisy);
    
    mm = length(xm);
ee = length(xd);
ss = length(x_source);
k=1;
for i = 1:dim:ss
    xs(k) = x_source(i);
    ys(k) = x_source(i+1);
    k = k +1;
end

ss = length(xs);

for i = 1:ee
    k = 1;
    for j = 1:ss
    r_mod1(i,j) = sqrt((xd(i)-xs(j))^2 + (yd(i)-ys(j))^2);
    V_jac(i,k) =  (1/2*pi)*log(r_mod1(i,j));
    V_jac(i,k+1) =  (S(j)/2*pi)*xs(j)/ (r_mod1(i,j))^2;
    V_jac(i,k+2) =  (S(j)/2*pi)*ys(j)/ (r_mod1(i,j))^2;
    k = k + 3;
    end
end

for i = 1:mm
    k = 1;
    for j = 1:ss
    r_mod2(i,j) = sqrt((xm(i)-xs(j))^2 + (ym(i)-ys(j))^2);
    Gs(i,k) =  -(1/2*pi)*1/(r_mod2(i,j));
    Gs(i,k+1) =  (S(j)/2*pi)*xs(j)/ (r_mod2(i,j))^3;
    Gs(i,k+2) =  (S(j)/2*pi)*ys(j)/ (r_mod2(i,j))^3;
    k = k + 3;
    end
end
[Ua S_a Va] = svd(A);
% A_inv = Va*inv(S_a)*Ua';
Sa_diag = diag(S_a);

%%%%%trial TSVD
% pp = 1;
% Sar = diag(Sa_diag(1:end-pp));
% Uar = Ua(:,1:end-pp);
% Var = Va(1:end-pp,:)';
% Ar = Uar*Sar*Var';
% Ar_inv = Var*inv(Sar)*Uar';


%%%%%%trial Tichnov
Dap = diag(Sa_diag./(Sa_diag.^2 + alph_tik));
Ar_inv = Va*Dap*Ua';
L = M*Ar_inv*P - N;
% L = M*rand(100,100)*P - N;
S_adjoint = V_jac + L*Gs;
% cond(Ar_inv)
% cond(S_adjoint)
[Us S_s Vs] = svd(S_adjoint);
SS_diag = diag(S_s);

legendInfo{kk} = ['rr=' num2str(rr) 'sym'];
% kk = kk + 1 ;
figure(10)
plot(log(SS_diag)); 
hold on;           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
n_source = length(S);
th2 = fan_range*rand(n_source,1);

    for i = 1:n_source
    R(i) = normrnd(rr,sigma);
    x_source1b(i) = R(i)*cos(th2(i));
    x_source2b(i) = R(i)*sin(th2(i));
    end
       k=1;
       for i = 1:dim:2*n_source-1
         x_source22(i) = x_source1b(k);
         x_source22(i+1) = x_source2b(k);
         k = k+1;
%     scatter( x_source(i), x_source(i+1)); hold on;
       end
     n_loc        = length(x_source22);

 %%%% Forward BEM method to solve boundary data phim 
    [ g2, pk2 ] = boundary_pk( r,theta,nx,ny,x,y,x_source22,S,sanity,dim );

    if(reg == 0)
    [ ff2,F1,A2,P2 ] = NtD( g2,x_source22,S,x,y,xx,yy,...
                                   nx,ny,xm,ym,lx,...
                                   ly,l,r,theta,sanity,dim,alph_tik);

    elseif (reg == 1)
        Ka = tril(ones(n,n),-1);
        K1 = ones(n,1);
        [ ff2,F1,A2,P2 ] = NtD_miu2( g2,x_source22,S,x,y,xx,yy,...
                                   nx,ny,xm,ym,lx,...
                                   ly,l,r,theta,sanity,dim,reg_miu,Ka,K1,...
                                   gmres_restart,gmres_tol,gmres_maxit);
    elseif (reg == 2)
        Ka = tril(ones(n,n),-1);
        K1 = ones(n,1);
        [ ff2,F1,A2,P2 ] = NtD_miu( g2,x_source22,S,x,y,xx,yy,...
                                   nx,ny,xm,ym,lx,...
                                   ly,l,r,theta,sanity,dim,reg_miu,Ka,K1,...
                                   gmres_restart,gmres_tol,gmres_maxit);
                                    
    end

%     for m=1:length(xm)                               
    for m=1:length(xd)
          F2 = integrals(l,nx,ny,xx,yy,xd(m),yd(m),lx,ly);
          f_data(m) =  2.0*(F2(2,:)*ff' - F2(1,:)*g');
          M2(m,:) = F2(2,:);
          N2(m,:) = F2(1,:);
    end     
    if (full_partial == 1)
    for m=1:length(xm)                               
%     for m=1:length(xd)
          F22 = integrals(l,nx,ny,xx,yy,xm(m),ym(m),lx,ly);
          f_data2(m) =  2*(F22(2,:)*g' - F22(1,:)*ff2');
    end
    end
    
    f_measure_noisy2 = awgn(f_data,snr);
    norm(f_measure_noisy2,2)
    
    figure(1)
    plot(x,y); hold on;
    for i = 1:dim:2*n_source
        scatter(x_source22(i),x_source22(i+1));
        hold on;
    end
    figure(2)
    plot(f_data); hold on; plot(f_measure_noisy);
    
    mm = length(xm);
ee = length(xd);
ss = length(x_source22);
k=1;
for i = 1:dim:ss
    xs(k) = x_source22(i);
    ys(k) = x_source22(i+1);
    k = k +1;
end

ss = length(xs);

for i = 1:ee
    k = 1;
    for j = 1:ss
    r_mod3(i,j) = sqrt((xd(i)-xs(j))^2 + (yd(i)-ys(j))^2);
    V_jac2(i,k) =  (1/2*pi)*log(r_mod3(i,j));
    V_jac2(i,k+1) =  (S(j)/2*pi)*xs(j)/ (r_mod3(i,j))^2;
    V_jac2(i,k+2) =  (S(j)/2*pi)*ys(j)/ (r_mod3(i,j))^2;
    k = k + 3;
    end
end

for i = 1:mm
    k = 1;
    for j = 1:ss
    r_mod4(i,j) = sqrt((xm(i)-xs(j))^2 + (ym(i)-ys(j))^2);
    Gs2(i,k) =  -(1/2*pi)*1/(r_mod4(i,j));
    Gs2(i,k+1) =  (S(j)/2*pi)*xs(j)/ (r_mod4(i,j))^3;
    Gs2(i,k+2) =  (S(j)/2*pi)*ys(j)/ (r_mod4(i,j))^3;
    k = k + 3;
    end
end
[Ua2 S_a2 Va2] = svd(A2);
% A_inv = Va*inv(S_a)*Ua';
Sa_diag2 = diag(S_a2);

%%%%%trial TSVD
% pp = 1;
% Sar = diag(Sa_diag(1:end-pp));
% Uar = Ua(:,1:end-pp);
% Var = Va(1:end-pp,:)';
% Ar = Uar*Sar*Var';
% Ar_inv = Var*inv(Sar)*Uar';


%%%%%%trial Tichnov
Dap2 = diag(Sa_diag2./(Sa_diag2.^2 + alph_tik));
Ar_inv2 = Va2*Dap2*Ua2';
L2 = M2*Ar_inv2*P2 - N2;
% L = M*rand(100,100)*P - N;
S_adjoint2 = V_jac2 + L2*Gs2;
% cond(Ar_inv)
% cond(S_adjoint)
[Us2 S_s2 Vs2] = svd(S_adjoint2);
SS_diag2 = diag(S_s2);

legendInfo{kk+1} = ['rr=' num2str(rr) 'asym'];
kk = kk + 2 ;
figure(10)
plot(log(SS_diag2),'*'); 
hold on;          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
end
  legend(legendInfo);
 

  
  

% for rr = 0.85;         % radius range with normal distribution
% 
% 
% 
% 
% end   
  
% title('4 pts, r=0.85, Tichnov, Symmetric ')
% ylim([-60 20])
% save('jac_Symmetric_r085.mat','SS_diag')


    % Get Data Variance Covariance matrix
%     Sig = (1/length(f_measure_noisy))*f_measure_noisy'*f_measure_noisy;
%     EigSig = eig(Sig);
    %%%% End of solving boundary data using BEM  
    
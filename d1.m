clc;
clear all;
close all;


global alph snr  
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

%%control pamameters for gmres solver:
gmres_restart = 25;
gmres_tol = 10^(-10); 
gmres_maxit = 1000;
%%%% Use geo_generate routine to generate artificial data
%%%% Load geometry 


load('D:\1 Research matlab\geo.mat');
    dim = 2;
   k_count = 1;
    for k_angle=0.001:0.01:0.2
  reg = 0;            % 0 for Tichnov NtD, 1 for miu regularised Ntd
                    % 2 for squared miu with kernel
reg_miu =0;      
    S = [1 -1];
    n_source = length(S);
    ang1 = -pi; small_ang = k_angle;
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


% Control parameters:
sanity = 0;        %% Dummy, change it
type = 'N';     % 1.Boundary type switch,  D--> Dirichlet, N --> Neumann.  
constrt = 1;
method = 1;    % 1 for fmincon, 2 for simulated annealing, 3 for hybrid
lambda_inv = 0.0; % Tichnov Regularisation for the inverse problem
weight1 = 1;  %weight for constraint SA, weight1 for equality
weight2 = 1;  %weight for constraint SA, weight2 for inequality
%%%% Load geometry 
load('geo.mat')

%%%% Load data
    load('data.mat');
    n_tot = n_source ;

%%%% Initial Guess for the minimization solver:
    x0 = rand(3*n_source,1);
%%%%%%%%%%%%%%%%%%%
f_c = @(x) circle_con(x,n_source);    
f = @(source_sol)BEM_getR_Neumann( source_sol,n_source,n_loc,x,y,xx,yy,nx,...
                                ny,xm,ym,lx,ly,l,r,theta,sanity,xd,yd,...
                                f_measure_noisy,dim,lambda_inv,alph_tik);

                          

%%%%%%%%%% Set up constraints for optimization %%%%%%%%%%

%%% 1. nutral source strength constraint
Aa = [(ones(n_source,1))' (zeros(2*n_source,1))'];    
ba = [0];

%%% 2. upper and lower bound for source strength constraint
K=2*max(S);
lb = [-K*ones(n_source,1);[]]';
ub = -lb;


%%%%%%%%%% Min solver %%%%%%%%%%%
opts = optimoptions(@fminunc,'PlotFcns',{@optimplotfval,@optimplotx},...
    'MaxIterations',90000,'MaxFunctionEvaluations',90000,...
    'OptimalityTolerance',1e-16,'StepTolerance',1e-16);
problem = createOptimProblem('fmincon','x0',x0, ...
    'objective',f,'Aeq',Aa,'beq',ba,'lb',lb,'ub',ub,'nonlcon',f_c,...
    'options',opts);


% Gradient methods:
[x_invsoln,fval,exitflag,output] = fmincon(problem);


S_inv = x_invsoln(1:n_source);
loc_inv = x_invsoln(1+n_source:end);
save('inv_soln.mat','loc_inv','S_inv')
x_plotinv = loc_inv(1:dim:end-1);
y_plotinv = loc_inv(2:dim:end);


for i = 1:n_source    
rms_xy(i,k_count) =((x_plot(i) - x_plotinv(i))^2 + (y_plot(i) - y_plotinv(i))^2);
rms_S(i,k_count) = ((1 - abs(S_inv(i)))^2);
end
rms_S = sqrt((1/n_source)*sum(rms_S));
rms_xy = sqrt((1/n_source)*sum(rms_xy));

k_count = k_count +1;

    end
k_angle=0.001:0.01:0.2;
save('trial1','k_angle','rms_S','rms_xy');
figure(5)
plot(k_angle,log(rms_S),'*');
xlabel('distance between the two sources')
ylabel('log round mean square error of intensities')
figure(6)
plot(k_angle,log(rms_xy),'*');
xlabel('distance between the two sources')
ylabel('log round mean square error of locations')







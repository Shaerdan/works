clc;
clear all;
close all;

% Control parameters:
global alph snr 
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

if (type =='N' )
f_c = @(x) circle_con(x,n_source);    
f = @(source_sol)BEM_getR_Neumann( source_sol,n_source,n_loc,x,y,xx,yy,nx,...
                                ny,xm,ym,lx,ly,l,r,theta,sanity,xd,yd,...
                                f_measure_noisy,dim,lambda_inv,alph_tik);

end                           
if(type =='D')
 load('NtD.mat') ;
 phim = g;       
 f = @(S,x_plot,y_plot)BEM_getR_Dirichlet(phim,S,x,n_source,n_loc,...
                            x,y,xx,yy,nx,ny,xm,ym,lx,...
                            ly,l,r,theta,sanity,xd,yd,...
                            phim_measure_noisy,dim,lambda1,lambda2);       

end

%%%%%%%%%% Set up constraints for optimization %%%%%%%%%%

%%% 1. nutral source strength constraint
Aa = [(ones(n_source,1))' (zeros(2*n_source,1))'];    
ba = [0];

%%% 2. upper and lower bound for source strength constraint
K=2*max(S);
lb = [-K*ones(n_source,1);[]]';
ub = -lb;


%%%%%%%%%% Min solver %%%%%%%%%%%
if (method == 1)
opts = optimoptions(@fminunc,'PlotFcns',{@optimplotfval,@optimplotx},...
    'MaxIterations',90000,'MaxFunctionEvaluations',90000,...
    'OptimalityTolerance',1e-16,'StepTolerance',1e-16);
if (constrt == 1)
problem = createOptimProblem('fmincon','x0',x0, ...
    'objective',f,'Aeq',Aa,'beq',ba,'lb',lb,'ub',ub,'nonlcon',f_c,...
    'options',opts);
else
    problem = createOptimProblem('fmincon','x0',x0, ...
    'objective',f,...
    'options',opts);
end
% Gradient methods:
[x_invsoln,fval,exitflag,output] = fmincon(problem);
end

if (method == 2)
options = saoptimset('PlotFcns',{@saplotbestx,...
          @saplotbestf,@saplotx,@saplotf});
% Simulated Annealing:
[x_invsoln,fval,exitflag,output] = simulannealbnd(f,x0,lb,ub,options);
end

if (method ==3 )
% Hybrid method:    
opts = optimoptions(@fminunc,'PlotFcns',{@optimplotfval,@optimplotx},...
    'MaxIterations',60000,'MaxFunctionEvaluations',60000,...
    'OptimalityTolerance',1e-16,'StepTolerance',1e-16);
if (constrt == 1)
problem = createOptimProblem('fmincon','x0',x0, ...
    'objective',f,'Aeq',Aa,'beq',ba,'lb',lb,'ub',ub,'nonlcon',f_c,...
    'options',opts);
else
    problem = createOptimProblem('fmincon','x0',x0, ...
    'objective',f,...
    'options',opts);
end
% First run gradient methods:
[x_invsoln,fval,exitflag,output] = fmincon(problem);
% Simulated Annealing: 
x0 =x_invsoln;
options = saoptimset('PlotFcns',{@saplotbestx,...
          @saplotbestf,@saplotx,@saplotf});
[x_invsoln,fval,exitflag,output] = simulannealbnd(f,x0,lb,ub,options);
end

% if (method == 4)
%         x01 = rand(n_tot,1);
%         x02 = -x01;
%         x03 = rand(n_tot,1);
%         x04 = -x03;
% 
%      opts = optimoptions(@fminunc,'PlotFcns',{@optimplotfval,@optimplotx},...
%     'MaxIterations',10000,'MaxFunctionEvaluations',10000,...
%     'OptimalityTolerance',1e-18,'StepTolerance',1e-18);
% 
%     problem01 = createOptimProblem('fmincon','x0',x01, ...
%     'objective',f,'Aeq',Aa,'beq',ba,'lb',lb,'ub',ub,'nonlcon',f_c,...
%     'options',opts);
%     problem02 = createOptimProblem('fmincon','x0',x02, ...
%     'objective',f,'Aeq',Aa,'beq',ba,'lb',lb,'ub',ub,'nonlcon',f_c,...
%     'options',opts);
%     problem03 = createOptimProblem('fmincon','x0',x03, ...
%     'objective',f,'Aeq',Aa,'beq',ba,'lb',lb,'ub',ub,'nonlcon',f_c,...
%     'options',opts);
%     problem04 = createOptimProblem('fmincon','x0',x04, ...
%     'objective',f,'Aeq',Aa,'beq',ba,'lb',lb,'ub',ub,'nonlcon',f_c,...
%     'options',opts);
% 
% 
% [x_invsoln01,fval01,exitflag01,output01] = fmincon(problem01);
% figure(4)
% plotsoln( x_invsoln01,n_source,dim,x,y,x_source )
% 
% [x_invsoln02,fval02,exitflag02,output02] = fmincon(problem02);
% figure(5)
% plotsoln( x_invsoln02,n_source,dim,x,y,x_source )
% 
% [x_invsoln03,fval03,exitflag03,output03] = fmincon(problem03);
% figure(6)
% plotsoln( x_invsoln03,n_source,dim,x,y,x_source )
% 
% [x_invsoln04,fval04,exitflag04,output04] = fmincon(problem04);
% figure(7)
% plotsoln( x_invsoln04,n_source,dim,x,y,x_source )
% 
% end

if (method ~= 4 )
S_inv = x_invsoln(1:n_source);
loc_inv = x_invsoln(1+n_source:end);
save('inv_soln.mat','loc_inv','S_inv')
x_plotinv = loc_inv(1:dim:end-1);
y_plotinv = loc_inv(2:dim:end);

close all;
figure(1)
% subplot(2,1,1)
    plot(x,y); hold on;
    k = 1;
    for i = 1:n_source
        scatter(x_plot(i),y_plot(i));
        text1 = S(i);
        text1 = num2str(text1);
        text1 = cellstr(text1);
        text_fin =text(x_plot(i)+0.05,y_plot(i)+0.05,text1);
        fontsize = text_fin.FontSize;
        text_fin.FontSize = 4;
        hold on;
    end
    
figure(2)
plot(x,y);
hold on;
k=1;
for i = 1:n_source
scatter(x_plotinv(i),y_plotinv(i),'k','*')
        text1 = S_inv(i);
        text1 = num2str(text1);
        text1 = cellstr(text1);
        text_fin =text(x_plotinv(i)+0.05,y_plotinv(i)+0.05,text1);
        fontsize = text_fin.FontSize;
        text_fin.FontSize = 4;
hold on;
end


%%%%%%%% Returning the data %%%%%%%%

   [ g, gk ] = boundary_pk( r,theta,n_source,nx,ny,x,y,S_inv,x_plotinv,y_plotinv,sanity );
    [ ff,F1,A,P ] = NtD( g,xx,yy,...
                         nx,ny,xm,ym,lx,...
                         ly,l,alph_tik);
    for m=1:length(xd)
          F = integrals(l,nx,ny,xx,yy,xd(m),yd(m),lx,ly);
          f_data_reconst(m) =  2*(F(2,:)*ff' - F(1,:)*g');
          M(m,:) = F(2,:);
          N(m,:) = F(1,:);
    end  
       [ ve_reconst ] = get_ve( xd,yd,S_inv,x_plotinv,y_plotinv );
       f_data_reconst = f_data_reconst + ve_reconst;
       nrm = norm(f_data_reconst - f_measure_noisy,2)^2;
       figure(3)
           theta_sample=-pi:(2*pi/(m-1)):pi;
       plot(theta_sample,f_data_reconst,'*','DisplayName',...
    sprintf('%s', 'reconstructed')); hold on; plot(theta_sample,f_data,'DisplayName',...
    sprintf('%s', 'Actual'));
    legend('show');
%     title(['Noise Amplification Factor']); 



%%%%%% Rearranging Solution to Match the exact solution setting %%%%%%%
eps0 = 0.05;
eps1 = 0.5;
for i = 1:n_source
    for j = 1:n_source
    test = abs(x_plot(i) - x_plotinv(j));
    if (test < eps0)
        x_plotinv1(i) = x_plotinv(j);
    end
    end
end

for i = 1:n_source
    for j = 1:n_source
    test = abs(y_plot(i) - y_plotinv(j));
    if (test < eps0)
        y_plotinv1(i) = y_plotinv(j);
    end
    end
end

for i = 1:n_source
    for j = 1:n_source
    test = abs(S(i) - S_inv(j));
    if (test < eps1)
        S_inv1(i) = S_inv(j);
    end
    end
end

%%%%% Quantified Error for reconstruction %%%%%%



for i = 1:n_source    
rms_xy(i) =((x_plot(i) - x_plotinv1(i))^2 + (y_plot(i) - y_plotinv1(i))^2);
rms_S(i) = ((1 - abs(S_inv1(i)))^2);
end
rms_S = sqrt((1/n_source)*sum(rms_S));
rms_xy = sqrt((1/n_source)*sum(rms_xy));

Tc = clusterdata([x_plot; x_plotinv'],0.1);
end
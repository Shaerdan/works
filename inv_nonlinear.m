clc;
clear all;
close all;

% Control parameters:
global alph snr 
sanity = 0;        %% Dummy, change it
type = 'N';     % 1.Boundary type switch,  D--> Dirichlet, N --> Neumann.  
constrt = 1;
method = 1;    % 1 for fmincon, 2 for simulated annealing, 3 for hybrid
lambda_inv = 0.0001; % Tichnov Regularisation for the inverse problem
weight1 = 1;  %weight for constraint SA, weight1 for equality
weight2 = 1;  %weight for constraint SA, weight2 for inequality
%%%% Use geo_generate routine to generate artificial data
%%%% Load geometry 
load('geo.mat')

%%%% Use data_generate routine to generate artificial data
%%%% Load data
    load('data.mat');
    n_tot = 3*n_source;
%     x0 = zeros(n_tot,1);
%     n_source = 2;
%     n_loc = n_source*2;
%     n_tot = n_source + n_loc;
%%%% Initial Guess for the minimization solver:
    x0 = rand(n_tot,1);
%     n_source = length(x0)/3;
%     n_loc = n_source*2;

%         x0 = [0;0;0;0;0;0;0;0;0];
%     n_source = 3;
%     n_loc = 6;
%%%%% Finished

%%%% R(r) Minimization Routine: 
% [ R ] = BEM_getR( x_source,S,x,y,xx,yy,nx,ny,xm,ym,lx,...
%                                 ly,l,r,theta,sanity,xd,yd,...
%                                 phim_measure_noisy);

%%%%%%%%%%%%%%%%%%%

if (type =='N' && reg ==0 )
f_c = @(x_plot,y_plot) circle_con(x_plot,y_plot);    
f = @(S,x_plot,y_plot)BEM_getR_Neumann( S,x_plot,y_plot,n_source,n_loc,...
                                x,y,xx,yy,nx,...
                                ny,xm,ym,lx,ly,l,r,theta,sanity,xd,yd,...
                                f_measure_noisy,dim,lambda_inv,alph_tik);
elseif (type =='N' && reg ==1 )  
 f_c = @(x_plot,y_plot) circle_con(x_plot,y_plot);   
 f = @(S,x_plot,y_plot)BEM_getR_Neumann_miu( S,x_plot,y_plot,n_source,n_loc,x,y,xx,yy,...
       nx,ny,xm,ym,lx,ly,l,r,theta,sanity,xd,yd,...
       f_measure_noisy,dim,lambda_inv,reg_miu,K1,gmres_restart,gmres_tol,gmres_maxit);
end                           
if(type =='D')
 load('NtD.mat') ;
 phim = g;       
 f = @(S,x_plot,y_plot)BEM_getR_Dirichlet(phim,S,x_plot,y_plot,x,n_source,n_loc,...
                            x,y,xx,yy,nx,ny,xm,ym,lx,...
                            ly,l,r,theta,sanity,xd,yd,...
                            phim_measure_noisy,dim,lambda1,lambda2);       

end
% [x_invsoln,fval,exitflag,output] = fminsearch(f,x0);
% options = optimoptions(@fminunc,'PlotFcns',@optimplotfval,'Algorithm',...
%     'quasi-newton','OptimalityTolerance',1e-16);
% [x_invsoln,fval,exitflag,output] = fminunc(f,x0,options);
Aa = [(ones(n_source,1))' (zeros(2*n_source,1))'];    % Make it flexible for more than 2 sources
ba = [0];
K=2;
% lb = [-K*ones(n_source,1); -ones(2*n_source,1)]';
lb = [-K*ones(n_source,1);[]]';
ub = -lb;
%%%%%%%%% ADD INEQUALITY BOUND%%%%%%%
if (method == 1)
opts = optimoptions(@fminunc,'PlotFcns',{@optimplotfval,@optimplotx},...
    'MaxIterations',60000,'MaxFunctionEvaluations',60000,...
    'OptimalityTolerance',1e-30,'StepTolerance',1e-30);
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
[S_inv,x_plotinv,y_plotinv,fval,exitflag,output] = fmincon(problem);
% Simulated Annealing:
end
if (method == 2)
options = saoptimset('PlotFcns',{@saplotbestx,...
          @saplotbestf,@saplotx,@saplotf});
[S_inv,x_plot,y_plot,fval,exitflag,output] = simulannealbnd(f,x0,lb,ub,options);
end
if (method ==3 )
opts = optimoptions(@fminunc,'PlotFcns',{@optimplotfval,@optimplotx},...
    'MaxIterations',6000,'MaxFunctionEvaluations',60000,...
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
[S_inv,x_plot,y_plot,fval,exitflag,output] = fmincon(problem);
% Simulated Annealing: 
x0 =x_invsoln;
options = saoptimset('PlotFcns',{@saplotbestx,...
          @saplotbestf,@saplotx,@saplotf});
[x_invsoln,fval,exitflag,output] = simulannealbnd(f,x0,lb,ub,options);
end

if (method == 4)
        x01 = rand(n_tot,1);
        x02 = -x01;
        x03 = rand(n_tot,1);
        x04 = -x03;

     opts = optimoptions(@fminunc,'PlotFcns',{@optimplotfval,@optimplotx},...
    'MaxIterations',10000,'MaxFunctionEvaluations',10000,...
    'OptimalityTolerance',1e-18,'StepTolerance',1e-18);

    problem01 = createOptimProblem('fmincon','x0',x01, ...
    'objective',f,'Aeq',Aa,'beq',ba,'lb',lb,'ub',ub,'nonlcon',f_c,...
    'options',opts);
    problem02 = createOptimProblem('fmincon','x0',x02, ...
    'objective',f,'Aeq',Aa,'beq',ba,'lb',lb,'ub',ub,'nonlcon',f_c,...
    'options',opts);
    problem03 = createOptimProblem('fmincon','x0',x03, ...
    'objective',f,'Aeq',Aa,'beq',ba,'lb',lb,'ub',ub,'nonlcon',f_c,...
    'options',opts);
    problem04 = createOptimProblem('fmincon','x0',x04, ...
    'objective',f,'Aeq',Aa,'beq',ba,'lb',lb,'ub',ub,'nonlcon',f_c,...
    'options',opts);


[S,x_plot,y_plot,fval01,exitflag01,output01] = fmincon(problem01);
figure(4)
plotsoln( x_invsoln01,n_source,dim,x,y,x_source )

[S,x_plot,y_plot,fval02,exitflag02,output02] = fmincon(problem02);
figure(5)
plotsoln( x_invsoln02,n_source,dim,x,y,x_source )

[S,x_plot,y_plot,fval03,exitflag03,output03] = fmincon(problem03);
figure(6)
plotsoln( x_invsoln03,n_source,dim,x,y,x_source )

[S,x_plot,y_plot,fval04,exitflag04,output04] = fmincon(problem04);
figure(7)
plotsoln( x_invsoln04,n_source,dim,x,y,x_source )

end

if (method ~= 4 )
S_inv = x_invsoln(1:n_source);
loc_inv = x_invsoln(1+n_source:end);
save('inv_soln.mat','loc_inv','S_inv')

close all;
figure(1)
    subplot(4,2,1)
    plot(x,y); hold on;
    scatter(xd,yd,'*')
    k = 1;
    for i = 1:dim:2*n_source
        scatter(x_source(i),x_source(i+1));
        text1 = S(k);
        text1 = num2str(text1);
        text1 = cellstr(text1);
        k = k + 1;
        text_fin =text(x_source(i)+0.05,x_source(i+1)+0.05,text1);
        fontsize = text_fin.FontSize;
        text_fin.FontSize = 4;
        hold on;
    end
    subplot(4,2,2)
    plot(f_data); hold on; plot(f_measure_noisy);


subplot(4,2,3)
plot(x,y);
hold on;
k=1;
for i = 1:dim:length(loc_inv)
scatter(loc_inv(i),loc_inv(i+1),'k','*')
        text1 = S_inv(k);
        text1 = num2str(text1);
        text1 = cellstr(text1);
        k = k + 1;
        text_fin =text(x_source(i)+0.05,x_source(i+1)+0.05,text1);
        fontsize = text_fin.FontSize;
        text_fin.FontSize = 4;
hold on;
end

figure(2)
plot(S); hold on; plot(S_inv)

end
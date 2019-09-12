clc;
clear all;
close all;

% Control parameters:
sanity = 0;        %% Dummy, change it
type = 'N';     % 1.Boundary type switch,  D--> Dirichlet, N --> Neumann.  
constrt = 1;
method = 1;    % 1 for fmincon, 2 for simulated annealing, 3 for hybrid
solver = 'tichnov and linear'; % tichnov and linear or tichnov

alph_tik2 = 10^(-1); % Tichnov parameter for the final inverse 
alph_linear = 10^12; % Linear constraint paramater
alph_deflate = 10^12; % Deflation for NtD map;
%%%% Use geo_generate routine to generate artificial data
%%%% Load geometry 
load('geo.mat')
load('geo2.mat')

%%%% Use data_generate routine to generate artificial data
%%%% Load data
load('datafixloc.mat')

n_mesh = n1;
n_e = length(xd);
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%     S = source_sol(1:n_source);
%     x_source = source_sol(n_source+1:end);
        
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


%%%%%Tichnov Regularisation:
[Ua, S_a, Va] = svd(A);
Sa_diag = diag(S_a);
Dap = diag(Sa_diag./(Sa_diag.^2 + 10^(-5)));
Ar_inv = Va*Dap*Ua';
z = Ar_inv*b;
ff = z';
% test0 =inv(P)*A*ff';


%%%%% NtD map Gamma regularized %%%%%
% GAMMA_reg = A'*A  +...
%      alph_deflate*((inv(P)*A)'*eye(n_mesh,1)*eye(1,n_mesh)*(inv(P)*A));
% ff =  GAMMA_reg\(A'*b);
% ff = ff';
% sum(ff)
Gamma = Ar_inv*P;

[ Me,Ne ] = fixloc_getMeNe( n_e,xd,yd,nx,ny,l,xx,yy,lx,ly );

 G = (Me*Gamma - Ne);

 C = fixloc_getC( n_mesh,n_source,xm,ym,x_plot,y_plot);
 
g_C = C*S;
k = 1;
for i = 1:length(g_C)-1
    g_Cm(k) = -(g_C(i) + g_C(i+1))/2;
    k = k + 1 ;
end
 
 Se = fixloc_getSe(  n_e,n_source,xd,yd,x_plot,y_plot );

  
 Gain = (Se + 2*G*C);
 cond(Gain)
 
 [U_G, S_G, V_G ] = svd(Gain);
 figure(1)
 plot(diag(log(S_G)),'*'); hold on;
 
%  S_G = diag(S_G);
SG_diag = diag(S_G); 

nn0 = length(SG_diag);
nn1 = n_source - nn0;

if (strcmp(solver,'tichnov'))
%  s_Ginv = (SG_diag./(SG_diag.^2 + alph_tik2));
 %%% take the form of Wiener filter:
 f_wiener = SG_diag./(SG_diag.^2 + alph_tik2);
 
  plot(log(f_wiener),'+'); hold off;

%  coef1 = f_wiener./s_Ginv;
 ub_prod = U_G'*f_measure_noisy';
 
 for i = 1:length(f_wiener)
     S_inv(:,i)=f_wiener(i)*ub_prod(i)*V_G(:,i);
 end   
 S_inv = sum(S_inv,2);
elseif (strcmp(solver,'tichnov and linear'))
    P_linear = [ones(1,n_source);zeros(n_source-1,n_source)];
 Gain_reg = Gain'*Gain + alph_tik2*eye(n_source,n_source) +...
     alph_linear*(P_linear'*P_linear);
 S_inv = Gain_reg\(Gain'*f_measure_noisy');
  f_wiener = SG_diag./(SG_diag.^2 + alph_tik2);

 Gain_reg_sing = svd(Gain_reg);
    plot(log(Gain_reg_sing),'+'); hold on;
    plot(log(1./f_wiener),'o')

end

sum(S_inv)


 figure(3)
 plot(S_inv); hold on; plot(S)
 
%  Trial = Gain*V_G*s_Ginv'*U_G';
 f_reconstrct = Gain*S_inv;
 figure(4)
 plot(f_measure_noisy); hold on; plot(f_reconstrct,'o');
 figure(5)
 plot(ff)
%  
 figure(6)
trisurf(tri,x_plot,y_plot,S);
caxis([-1 1])
figure(7)
trisurf(tri,x_plot,y_plot,S_inv);
caxis([-1 1])

figure(8)
Forward = Gain'*Gain;
plot(log(svd(Forward)),'*'); hold on;
plot(log(svd(Gain_reg)),'o')

save('gain.mat','Gain')
 
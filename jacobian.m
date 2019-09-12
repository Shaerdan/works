clc;
clear all;
close all;
load('geo.mat');
load('data.mat');

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
    r_mod(i,j) = sqrt((xm(i)-xs(j))^2 + (ym(i)-ys(j))^2);
    Vs(i,k) =  (1/2*pi)*log(r_mod(i,j));
    Vs(i,k+1) =  (S(j)/2*pi)*xs(j)/ (r_mod(i,j))^2;
    Vs(i,k+2) =  (S(j)/2*pi)*ys(j)/ (r_mod(i,j))^2;
    k = k + 3;
    end
end

for i = 1:mm
    k = 1;
    for j = 1:ss
    r_mod(i,j) = sqrt((xm(i)-xs(j))^2 + (ym(i)-ys(j))^2);
    Gs(i,k) =  -(1/2*pi)*1/(r_mod(i,j));
    Gs(i,k+1) =  (S(j)/2*pi)*xs(j)/ (r_mod(i,j))^3;
    Gs(i,k+2) =  (S(j)/2*pi)*ys(j)/ (r_mod(i,j))^3;
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
S_adjoint = Vs + L*Gs;
cond(Ar_inv)
cond(S_adjoint)
[Us S_s Vs] = svd(S_adjoint);
SS_diag = diag(S_s);
figure(10)
plot(log(SS_diag)); hold on;
% title('4 pts, r=0.85, Tichnov, Symmetric ')
ylim([-60 20])
save('jac_Symmetric_r085.mat','SS_diag')
% figure(11)
% plot(log(diag(Dap)));hold on;plot(log(Sa_diag))
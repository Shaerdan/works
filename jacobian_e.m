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
    r_mod(i,j) = sqrt((xd(i)-xs(j))^2 + (yd(i)-ys(j))^2);
    Sigma(j) = S(j)/(r_mod(i,j))^2;
    end
    Ve(i,2*i-1) = (xd(i)/(2*pi))*sum(Sigma);
    Ve(i,2*i) = (yd(i)/(2*pi))*sum(Sigma);
end

for i = 1:ee
    k = 1;
    for j = 1:ss
    r_mod(i,j) = sqrt((xd(i)-xs(j))^2 + (yd(i)-ys(j))^2);
    Sigma(j) = S(j)/(r_mod(i,j))^3;
    end
    Ge(i,2*i-1) = (xd(i)/(2*pi))*sum(Sigma);
    Ge(i,2*i) = (yd(i)/(2*pi))*sum(Sigma);
end

[U S_A V] = svd(A);
A_inv = V*inv(S_A)*U';
L = M*A_inv*P - N;
S_adjoint = Ve + L*Ge;

% for i = 1:mm
%     k = 1;
%     for j = 1:ss
%     r_mod(i,j) = sqrt((xm(i)-xs(j))^2 + (ym(i)-ys(j))^2);
%     Ge(i,k) =  -(1/2*pi)*1/(r_mod(i,j));
%     Ge(i,k+1) =  (S(j)/2*pi)*xs(j)/ (r_mod(i,j))^3;
%     Ge(i,k+2) =  (S(j)/2*pi)*ys(j)/ (r_mod(i,j))^3;
%     k = k + 3;
%     end
% end
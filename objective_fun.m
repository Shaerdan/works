function [J] = objective_fun(X)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
global n_fft an2 bn2 n_source theta n_sample 
global eps0 data_measure lam_J reg_switch eps eps_ball
k=1;
for i = 1:n_source
S(i) = X(k);
rs(i) = X(k+1);
ths(i) = X(k+2);
k = k+3;
end

for i = 1:n_source 
phi(:,i) = (S(i)/(2*pi))*log(sqrt(1+rs(i)^2-2*rs(i)*cos(theta-ths(i))));
end

phi = 2*sum(phi,2);

% y1 = fft(phi,n_sample);                             % Compute DFT of x
% m1 = abs(y1);                               % Magnitude
% p1 = unwrap(angle(y1));                     % Phase
% an1 = (2/(n_fft-1))*real(y1);
% bn1 = -(2/(n_fft-1))*imag(y1);

% J = norm(an1-an2,2)^2 + norm(bn1-bn2,2)^2;
if (reg_switch == 0)
J = norm(phi - data_measure,2)^2;
end


%%%%%% Regularizor %%%%%%
% if (reg_switch == 1)
% % bc integral grad phi^2:
% for i = 1:n_source
% fun = @(xx) (rs(i)*sin(xx-ths(i))./(1+rs(i)^2-2*rs(i)*cos(xx - ths(i)))).^2;
% reglr(i) = (1/2*pi)*integral(fun,-pi,pi);
% end
% reglr = (sum(reglr));
% 
% J = norm(phi - data_measure,2)^2;
% 
% J = J + lam_J*reglr;
% end

if (reg_switch == 1)
% bc integral grad phi^2:
for i = 1:n_source
fun = @(xx) (rs(i)*sin(xx-ths(i))./(1+rs(i)^2-2*rs(i)*cos(xx - ths(i)))).^2;
reglr(i) = (1/2*pi)*integral(fun,-pi,pi);
end
reglr = (sum(reglr));

J = norm(phi - data_measure,2)^2;

J = J + lam_J*reglr;
end

if (reg_switch == 9)
% full hessian regularisor
for i = 1:n_source
fun1 = @(xx) (rs(i)*sin(xx-ths(i))./(1+rs(i)^2-2*rs(i)*cos(xx - ths(i)))).^2;
fun2 = @(xx)((-8*rs(i)*cos(xx-ths(i))+2*(cos(xx-ths(i))).^2 + 1+rs(i)^2 )./(1+rs(i)^2-2*rs(i)*cos(xx - ths(i))).^2).^2;
fun3 = @(xx) ((2*rs(i)*((rs(i)^2 + 1)*cos(xx-ths(i)) - 2*rs(i)*(sin(xx-ths(i))).^2 - 2*rs(i)*(cos(xx-ths(i))).^2))./(-2*rs(i)*cos(xx-ths(i)) + rs(i)^2 + 1).^2).^2;
fun4 = @(xx) ((2*rs(i) - cos(xx-ths(i)))./(1+rs(i)^2-2*rs(i)*cos(xx - ths(i)))).^2;
end
reglr1 = (1/2*pi)*integral(fun1,-pi,pi);
reglr2 = (1/2*pi)*integral(fun2,-pi,pi);
reglr3 = (1/2*pi)*integral(fun3,-pi,pi);
reglr4 = (1/2*pi)*integral(fun3,-pi,pi);
reglr = reglr1 ;


J = norm(phi - data_measure,2)^2;

J = J + lam_J*reglr;
end


%%%%%% L1 norm bc regularizor (TV):
if (reg_switch == 2)

for i = 1:n_source
fun = @(xx) abs(log(1+rs(i)^2 - 2*rs(i)*cos(xx-ths(i))));
reglr(i) = (1/2*pi)*integral(fun,-pi,pi);
end
reglr = sum(reglr);

J = norm(phi - data_measure,2)^2;

J = J + lam_J*reglr;
end

if (reg_switch == 3)
for i = 1:n_source
% factor1 = (S(i)^2);    
factor1 = 1;    
l1 = @(xx,yy) (yy.^2+rs(i)^2-2*rs(i)*yy.*cos(xx - ths(i)))+eps0;
l2 =  @(xx,yy) (yy.^2+rs(i)^(-2) - 2*(yy/(rs(i)+eps0)).*cos(xx - ths(i)) )+eps0;
dfdr = @(xx,yy) (2*yy-2*rs(i)*cos(xx - ths(i)))./l1(xx,yy) + (2*yy-(2/(rs(i)+eps0))*cos(xx - ths(i)))./l2(xx,yy);
dfdth = @(xx,yy) (2*rs(i)*yy.*sin(xx - ths(i)))./l1(xx,yy) +(2*(yy/(rs(i)+eps0)).*sin(xx - ths(i)))./l2(xx,yy);
fun = @(xx,yy) yy.*(dfdr(xx,yy).^2 + ((1./(yy+eps0)).*dfdth(xx,yy)).^2);
% fun = @(xx,yy) f(xx,yy)./(1+eps*f(xx,yy));
% reglr(i) = factor1*(integral2(fun,-pi,ths(i),0,1) + integral2(fun,ths(i),pi,0,1));
% reglr(i) = factor1*(integral2(fun,-pi,pi,0,rs(i)-eps_ball) + integral2(fun,-pi,pi,rs(i)+eps_ball,1));
reglr(i) = factor1*(integral2(fun,-pi,ths(i)-eps_ball,0,1) + integral2(fun,ths(i)+eps_ball,pi,0,1));
end
reglr = sum(reglr);

J = norm(phi - data_measure,2)^2;

J = J + lam_J*reglr;
end
% figure(3)
% plot(theta,phi,theta,data_measure)
if (reg_switch == 4)
for i = 1:n_source
% factor1 = abs(S(i));   
factor1 = 1;    
l1 = @(xx,yy) (yy.^2+rs(i)^2-2*rs(i)*yy.*cos(xx - ths(i)))+eps0;
l2 =  @(xx,yy) (yy.^2+rs(i)^(-2) - 2*(yy/(rs(i))).*cos(xx - ths(i)))+eps0;
dfdr = @(xx,yy) (2*yy-2*rs(i)*cos(xx - ths(i)))./l1(xx,yy) + (2*yy-(2/(rs(i)))*cos(xx - ths(i)))./l2(xx,yy);
dfdth = @(xx,yy) (2*rs(i)*yy.*sin(xx - ths(i)))./l1(xx,yy) +(2*(yy/(rs(i))).*sin(xx - ths(i)))./l2(xx,yy);
fun = @(xx,yy) yy.*sqrt(dfdr(xx,yy).^2 + ((1./(yy)).*dfdth(xx,yy)).^2);
% fun = @(xx,yy) f(xx,yy)./(1+eps*f(xx,yy));
reglr(i) = factor1*(integral2(fun,-pi,pi,0,rs(i)-eps_ball) + integral2(fun,-pi,pi,rs(i)+eps_ball,1));
% reglr(i) = factor1*(integral2(fun,-pi,pi,0,rs(i)) + integral2(fun,-pi,pi,rs(i),1));
% reglr(i) = factor1*(integral2(fun,-pi,ths(i)-eps_ball,0,1) + integral2(fun,ths(i)+eps_ball,pi,0,1));
end
reglr = sum(reglr);
reglr = reglr;

J = norm(phi - data_measure,2)^2;

J = J + lam_J*reglr;
end

if (reg_switch == 5)
for i = 1:n_source
factor1 = 1;    
l1 = @(xx,yy) (yy.^2+rs(i)^2-2*rs(i)*yy.*cos(xx - ths(i)))+eps0;
l2 =  @(xx,yy) (yy.^2+rs(i)^(-2) - 2*(yy/(rs(i)+eps0)).*cos(xx - ths(i)))+eps0;
dfdr = @(xx,yy) (2*yy-2*rs(i)*cos(xx - ths(i)))./l1(xx,yy) + (2*yy-(2/(rs(i)+eps0))*cos(xx - ths(i)))./l2(xx,yy);
dfdth = @(xx,yy) (2*rs(i)*yy.*sin(xx - ths(i)))./l1(xx,yy) +(2*(yy/(rs(i)+eps0)).*sin(xx - ths(i)))./l2(xx,yy);
f = @(xx,yy) yy.*sqrt(dfdr(xx,yy).^2 + ((1./(yy+eps0)).*dfdth(xx,yy)).^2);
fun = @(xx,yy) f(xx,yy)./(1+eps*f(xx,yy));
reglr(i) = factor1*(integral2(fun,-pi,pi,0,1));
end
reglr = sum(reglr);

J = norm(phi - data_measure,2)^2;

J = J + lam_J*reglr;
end


if (reg_switch == 6)
    for i = 1:n_source 
    l1 = @(xx) (1/(2*pi))*log(sqrt(1+rs_real(i)^2 ...
    -2*rs_real(i)*cos(xx-ths_real(i))));
    end

end


end


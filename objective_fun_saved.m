function [J] = objective_fun(X)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
global n_fft an2 bn2 n_source theta n_sample data_measure lam_J reg_switch
global S rs ths
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

y1 = fft(phi,n_sample);                             % Compute DFT of x
m1 = abs(y1);                               % Magnitude
p1 = unwrap(angle(y1));                     % Phase
an1 = (2/(n_fft-1))*real(y1);
bn1 = -(2/(n_fft-1))*imag(y1);

% J = norm(an1-an2,2)^2 + norm(bn1-bn2,2)^2;
if (reg_switch == 0)
J = norm(phi - data_measure,2)^2;
end


%%%%%% Regularizor %%%%%%
if (reg_switch == 1)
% integral grad phi^2:
for i = 1:n_source
front_fun = ((1/2*pi)*S(i)*rs(i));
fun = @(xx) (sin(xx-ths(i))./(1+rs(i)^2-2*rs(i)*cos(xx - ths(i)))).^2;
reglr(i) = front_fun*integral(fun,-pi,pi);
end
reglr = sum(reglr)

J = norm(phi - data_measure,2)^2;

J = J + lam_J*reglr;
end

%%%%%% L1 norm regularizor (TV):
if (reg_switch == 2)

for i = 1:n_source
fun = @(xx) abs(log(1+rs(i)^2 - 2*rs(i)*cos(xx-ths(i))));
reglr(i) = (1/2*pi)*S(i)*integral(fun,-pi,pi);
end
reglr = sum(reglr)

J = norm(phi - data_measure,2)^2;

J = J + lam_J*reglr;
end


%%%%%% L1 norm of gradient regularizor:
if (reg_switch == 3)
% integral grad phi^2:
for i = 1:n_source
fun = @(xx) (sin(xx-ths(i))./(1+rs(i)^2-2*rs(i)*cos(xx - ths(i)))).^2;
reglr(i) = (1/2*pi)*S(i)*rs(i)*integral(fun,-pi,pi);
end
reglr = sum(reglr)

J = norm(phi - data_measure,2)^2;

J = J + lam_J*reglr;
end

if (reg_switch == 4)
% integral grad phi^2:
fun = @(x) tv_phi(x);
reglr =integral(fun,-pi,pi);

J = norm(phi - data_measure,2)^2;

J = J + lam_J*reglr^2;
end


% figure(3)
% plot(theta,phi,theta,data_measure)
end


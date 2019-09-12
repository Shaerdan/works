clc;
clear all;
close all;


global S rs ths k_source eps


S = [1 -1];
n_source = length(S);
ang_stp = (2*pi)/n_source;
% ths = 0;
ths = [0 0.0001];

k = 1;  
h_res = 0.1;
Rs =   0:h_res:0.99;  
for rs = 0:h_res:0.99



samples = 100;
level = (10^-1);

% for isam = 1:6
% clear phi_adjrs phi_adjS phi_adjths Adj;
n_sample = samples;
h_integral = 2*pi/n_sample;
kk = 1;

f1 = @(x,r,theta) -(1/(2*pi))*integrand_fourier_dS(x,r,theta);

f2 =@(x,r) integrand2( x,r );

eps = 10^(-6);

theta = -pi+eps:h_integral:pi-eps;
ll_theta = length(theta);

for k_source = 1:n_source
[ phi ] = analytical_dirichlet( theta,f1,eps );
end

y = fft(phi);                             % Compute DFT of x
m = abs(y);                               % Magnitude
p = unwrap(angle(y));                     % Phase
an = (2/99)*real(y);
a1(k) = an(2);
a2(k) = an(3);
bn = imag(y);



ff = (0:length(y)-1)*100/length(y);        % Frequency vector

% figure(1)
% subplot(2,1,1)
% plot(ff,m)
% title('Magnitude')
% ax = gca;
% ax.XTick = [15 40 60 85 100];
% 
% subplot(2,1,2)
% plot(ff,p*180/pi)
% title('Phase')
% ax = gca;
% ax.XTick = [15 40 60 85 100];



k = k + 1;

end

figure(2)
plot(Rs,abs(a1)); hold on; 
figure(3)
plot(Rs,abs(a2));
figure(4)
plot(Rs,abs(a1./a2));

ratio12 = abs(a1./a2);
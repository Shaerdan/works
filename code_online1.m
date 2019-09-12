% Adaptive Linear Predictive Coding
% This M file recursively calculated the coeficients of
% an unknown system. The value mu specifies the rate
% of convergence of the adaptive algorithm. mu must
% be less that the smallest eigenvalue of the unknown
% transfer function for the adaptive filter to be properly
% conditioned. len is the order of the filter.
% lin is the number of iterations to take place.
% x is the input signal to the system.
% n0, n1 are two independent noise sources.
% Author: Thomas Drumright
% Date: 6/11/98
clc;
clear all;
close all;
mu = .00000001; % Convergence factor
 % Filter length
a = 0;
mse = 0;
% Create input of the modeled system.
s = 0;
% for n=1:(lin);
%  s1(n) = 1*sin(2*pi.*n/abs(50*sin(pi*(2500.+n)/10000))).*sin(2*pi*n/200);
% end
M = csvread('C:\Users\SerdanSattar\Desktop\data\sp500.csv');
k=1;
for i = 17000:1:17200
ind_end = i;
data_x= M(1:ind_end,2);
s = data_x(1:end-1);
len = 5;
lin = ind_end-1; % Number of iterations
range = 100;
% Apply the adaptive algorithm
x = 0*(1:len)'; % Initial input
w = 0*(1:len); % Initial filter coefs.
% Calculate Nth iteration
for j = lin-range:lin
 y = w*x;
 a(j) = y;
 for n = len:-1:2
 x(n) = x(n-1);
 end
 x(1) = s(j);
 e = x(1)-y;
 mse(j+1) = e;
 w = w + mu*e*x';
end

y_predict(k) = w*x;
y_exact(k) = data_x(end);
y_prev(k) = data_x(end-1);
k = k+1;
end

figure(10)
plot(y_exact,'-*');hold on; plot(y_predict,'-o')

figure(11)
plot(y_prev); hold on; plot(y_predict)

figure(12)
plot(y_prev); hold on; plot(y_exact)


% Display the time domain of the filtered and
% unfiltered system.
% f = 0:length(a)-1;
% figure(1)
% hold on
% title('y(t)')
% 
% plot(f,a)
% l = 0:len-1;
% figure (2)
% hold on
% title('w(t)')
% plot (l,w)
% % Calculate and plot the mean squared error for
% % each iteration.
% num = 0:length(mse)-1;
% figure (3)
% hold on
% title('e(t)')
% plot(num,mse)
% % Calculate and plot the frequency spectrums of the
% % filtered and unfiltered signals. 
% W = fft(s,512);
% figure (4)
% hold on
% title('S(w)')
% t = 0:length(W)-1;
% plot(t,abs(W))
% figure (5)
% hold on
% title('s(t)')
% t = 0:length(s)-1;
% plot(t,s); W = fft(a,512);
% figure (6)
% hold on
% title('Y(w)')
% t = 0:length(W)-1;
% plot(t,abs(W)) ;W = fft(mse,512);
% figure (7)
% hold on
% title('E(w)')
% t = 0:length(W)-1;
% plot(t,abs(W))


% Display final coeficients and MSE
%disp('Estimated coeficients w =')
%disp(w')
%disp('Final error e =')
%disp(abs(e)) 
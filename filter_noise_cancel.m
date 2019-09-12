
clc;
clear all;
close all;

M = csvread('C:\Users\SerdanSattar\Desktop\data\sp500.csv');
x= M(:,2);


signal = x;

noise = randn(size(x));
nfilt = fir1(11,0.05); % Eleventh order lowpass filter
fnoise = filter(nfilt,1,noise); % Correlated noise data
d = x;


coeffs = nfilt.' -0.01; % Set the filter initial conditions.
mu = 0.0000005;          % Set the step size for algorithm updating.


lms = dsp.LMSFilter(12,'Method','Sign-Data LMS',...
   'StepSize',mu,'InitialConditions',coeffs);
[~,e] = step(lms,noise,d);
L = 500;
plot(0:L-1,signal(1:L),0:L-1,e(1:L));
title('Noise Cancellation by the Sign-Data Algorithm');
legend('Actual Signal','Result of Noise Cancellation',...
       'Location','NorthEast');
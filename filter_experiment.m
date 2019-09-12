clc;
clear all;
close all;
% openExample('dsp/PredictMeanSquareErrorForLMSFilterExample')


% step(obj,x);
% release(obj)

lms = dsp.LMSFilter;

num = fir1(31,0.5);
fir = dsp.FIRFilter('Numerator',num);  
iir = dsp.IIRFilter('Numerator',sqrt(0.75),...
        'Denominator',[1 -0.5]);
    M = csvread('C:\Users\SerdanSattar\Desktop\data\sp500.csv');
h1=1;
for i = 1 : 1000
h2=200+i;
    range = h1:h2;
    x = M(range,2);

% n = 0.1*randn(size(x));           
% d = x;
% 
% 
% l = 32;
% mu = 0.008;
% m  = 5;
% 
% lms = dsp.LMSFilter('Length',l,'StepSize',mu);
% [mmse,emse,meanW,mse,traceK] = msepred(lms,x,d,m);
% [simmse,meanWsim,Wsim,traceKsim] = msesim(lms,x,d,m);
% 
% nn = m:m:size(x,1);
% semilogy(nn,simmse,[0 size(x,1)],[(emse+mmse)...
% (emse+mmse)],nn,mse,[0 size(x,1)],[mmse mmse])
% title('Mean Squared Error Performance')
% axis([0 size(x,1) 0.001 10])
% legend('MSE (Sim.)','Final MSE','MSE','Min. MSE')
% xlabel('Time Index')
% ylabel('Squared Error Value')

order = 12;
b = fir1(order,0.8,'low');
d = filter(b,1,x(1:end-1));
mu = 0.9;
lms = dsp.LMSFilter(order+1,'StepSize',mu);
[y,e,w] = step(lms,x(1:end-1),d);

a = lpc(d,4);
est_x = filter([0 -a(2:end)],1,d);

error_est(i) = x(end)-est_x(end);

error_est2(i) = x(end)-x(end-1);
end

stem([b.' w]);

figure(4)
plot(error_est);

figure(5)
% e = x-est_x;
[acs,lags] = xcorr(e,'coeff');
plot(x(end-10:end)); hold on;
plot(est_x(end-10:end));

figure(6)
plot(x); hold on;
plot(est_x);
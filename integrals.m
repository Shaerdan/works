function [ F, delta ] = integrals( l,nx,ny,x,y,x0,y0,lx,ly )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

eps = 5*10^(-18);
A = l.^2;
B = 2*[-ny.*(x - x0) + nx.*(y - y0)].*l;
E = (x - x0).^2 + (y - y0).^2;
delta = 4*A.*E - B.^2;    

for i = 1: lx-1
        tpc = delta(i);
        tp1 = l(i)/(2*pi);
        tp2 = log(l(i));
        tp3 = B(i)/(2*A(i));
    if (abs(tpc) <= eps)
     F_1(i) = tp1 * (tp2 + (1+tp3)*log(abs(1+tp3)) ...
     - tp3 * log(abs(tp3)) - 1);
     F_2(i) = 0;
    else
     tea = E(i)/A(i);   
     tatan1 = atan((2*A(i) + B(i))/sqrt(tpc));
     tatan2 = atan( B(i)/sqrt(tpc));
     F_1(i) = (1/2) * tp1* ( 2*( tp2 -1 ) ... 
     - tp3*log(abs(tea)) + ...
     (1 + tp3)*log(abs(1 + 2*tp3 + tea))+((sqrt(tpc))/A(i))* ...
     (tatan1 - tatan2));
     F_2(i) = 2*tp1*(nx(i)*(x(i) - x0) + ny(i)*(y(i) - y0))...
         *(tatan1 - tatan2)/sqrt(tpc);
    end
     F = [F_1; F_2];
%      F_2(i) = l(i)*[nx(i)*(x(i)-x0)+ny(i)*(y(i)-y0)] 
end


function [J] = fourier_obj_fun(X)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
global k_star ak bk n_source lambda
k=1;
for i = 1:n_source
S(i) = X(k);
rs(i) = X(k+1);
ths(i) = X(k+2);
k = k+3;
end

co1 = 1/pi;
    for i = 1:k_star
        ak1(i) = -(S(1)/pi)*((rs(1)^i)/i)*cos(i*(ths(1)));
        ak2(i) = -(S(2)/pi)*((rs(2)^i)/i)*cos(i*(ths(2)));
        bk1(i) = -(S(1)/pi)*((rs(1)^i)/i)*sin(i*(ths(1)));
        bk2(i) = -(S(2)/pi)*((rs(2)^i)/i)*sin(i*(ths(2)));
        akh1(i) = -co1*((rs(1)^i)*sin(i*ths(1)));
        akh2(i) = -co1*((rs(2)^i)*sin(i*ths(2)));
        bkh1(i) = co1*((rs(1)^i)*cos(i*ths(1)));
        bkh2(i) = co1*((rs(2)^i)*cos(i*ths(2)));
    end
    
    aak = ak1+ak2;
    bbk = bk1+bk2;

J1 = (aak - ak).^2 + lambda*(akh1.^2 +akh2.^2);
J2 = (bbk - bk).^2 + lambda*(bkh1.^2 +bkh2.^2);


J = sum(J1) + sum(J2);

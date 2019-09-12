function [F] = auxilary_eqn(X)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
global k_star bk  lambda

I = X(1);
r = 0.7;
ths = X(3);

kk=1;
for i=1:k_star
rk(i) = r.^i;
rkm1(i) = r.^(i-1);
rkm2(i) = r.^(i-2);
r2km2(i) = r.^(2*i-2);
r2km1(i) = r.^(2*i-1);
r2k(i) = r.^(2*i);
sink(i) = sin(i*ths);
cosk(i) = cos(i*ths);
sin2k(i) = sin(2*i*ths);
cos2k(i) = cos(2*i*ths);
end
K = 1:k_star;



F(1) = I*sum(sin2k.*r2k./K) +pi*sum(cosk.*bk.*rk);
F(2) = I*sum(r2k./(K.^2)) - I*sum(cos2k.*r2k./(K.^2)) +...
    pi*sum(sink.*bk.*rk./K);
F(3) = I^2*sum(r2km1./K) - I^2*sum(cos2k.*r2km1./K) +...
    lambda*sum(K.*r2km1) + pi*I*sum(bk.*rkm1.*sink);


end


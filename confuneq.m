function [c,ceq] = confuneq(x)
% Nonlinear inequality constraints
global n_source
k=1;
for i = 1:n_source
    c(i) = -x(k+1);
    k = k + 3;
end
% Nonlinear equality constraints



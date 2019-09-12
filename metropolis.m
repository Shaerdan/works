function [ x_invsoln ] = metropolis( f,x0,h )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%   Metropolis algorithm to solve global minimum of function f,
%   x0 the initial guess

eps = 10^(-12);
n = length(x0);
f0 = f(x0);
step = -0.1;
h = 0.1;
grad = 0.001;
for T = 100:step:0.001
   T
   x1 = x0 - T*grad';
   for i = 1 : length(x1)
         grad(i) = (f(x1) - f(x0))/(x1(i) - x0(i));
   end
    dE = f(x1) - f0;
    if (dE <= eps)
        x0 = x1;
        f0 = f(x1);
        1
    else
    P = exp(-dE/T);
    r = rand;
       if ( P > r )
           x0 = x1;
           f0 = f(x1)
           0
       end
    end 

end
end


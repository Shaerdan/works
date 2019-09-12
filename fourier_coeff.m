function [an,bn] = fourier_coeff(phi)
global k_star
for i = 1:k_star
    fun1(xx) = cos(i*xx)*phi(xx);
    fun2(xx) = sin(i*xx)*phi(xx);
    ak(i) = integral(func1,-pi,pi);
    bk(i) = integral(func2,-pi,pi);
end


end
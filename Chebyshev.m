function [x,i] = Chebyshev(A,b,x0,M,tol,MaxEigs,MinEigs,iterNum)

d = (MaxEigs + MinEigs) / 2;
c = (MaxEigs - MinEigs) / 2;
x = x0;
r = b - A * x;
for i = 1:iterNum
    z = M*r;
    if (i == 1)
        p = z;
        alpha = 1/d;
    elseif (i == 2)
        beta = (1/2) * (c * alpha)^2;
        alpha = 1/(d - beta / alpha);
        p = z + beta * p;
    else
        beta = (c * alpha / 2)^2;
        alpha = 1/(d - beta / alpha);
        p = z + beta * p;
    end

    x = x + alpha * p;
    r = b - A * x;
    if (norm(r) < tol), break; end
end
end
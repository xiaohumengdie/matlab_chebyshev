function [x, m] = CG_Gear(A, b, x, M, res_tol, max_iter)

% CG method with different formular for alpha, by Chronogonlas Gear, 1992
% Correspond to Algoithm 2.3 in Gear's paper "s-step Iterative Methods
% for Symmetric Linear Systems"




s  = A*x;
r  = b - s;
z  = M*r;
rho_old = r'*z;
s  = z;
q  = A*s;

sigma = s'*q;
alpha  = rho_old / sigma;

x = x + alpha * s;
r = r - alpha * q;

for m = 1:max_iter 
    z = M*r;
    AZ= A * z;
    rho= r'*z;delta= AZ'*z;
    beta = rho/rho_old;
    sigma = delta-(beta*beta)*sigma;
    %sigma = delta-(beta*rho)/alpha;
    alpha = rho/sigma;rho_old = rho;
    s = z +beta*s;
    q = AZ+beta*q;
    x = x+alpha*s;
    r = r-alpha*q;
    if(norm(r,2)<res_tol)
        break
    end
end
end
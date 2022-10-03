%lanczos_based Algorithm
%A以及M都是对称负定
function [eigen, T] = lanczos_M(A, M)

[~,n]=size(A);
m = n;T=zeros(m,m);
q0=0;r0=ones(n,1);s0=M*r0;beta=0;
q1=r0/(sqrt(-r0'*s0));%更新q1
for j = 1:m
    p=M*q1;%更新p
    r=A*p-beta*q0;
    alpha=-r'*p;
    r=r-alpha*q1;%更新r
    s1=M*r;
    beta=sqrt(-r'*s1);
    q0=q1;q1=r/beta;
    %disp(beta)
    if(j==m)
    else
        T(j,j+1) = beta;
        T(j+1,j) = beta;
    end
    T(j,j)   = alpha;
end
eigen=eig(T);
end
format long ;
clc;clear;

tol=1e-6;cg2dpcoffdfac=0.51;

%A以及M_mitgcm,M_diagonal都是对称负定的
[A, b, x, M_mitgcm, M_mitgcm1, M_diagonal] = init_cg2d(cg2dpcoffdfac);

%M_ones为一个负的单位阵，对称负定
M_ones=-diag(ones(320,1));
[~,n]=size(A);
x0 = ones(n,1);x0 = x0/norm(x0,2);
iterNum = 1000000;

%也可以通过max(eig(A*M_mitgcm))以及min(eig(A*M_mitgcm))求取最大最小值
[eigen1, ~] = lanczos_M(A,M_ones);
[eigen2, ~] = lanczos_M(A,M_diagonal);%最大值可能不精确，导致结果发散。
[eigen3, ~] = lanczos_M(A,M_mitgcm);
[eigen4, ~] = lanczos_M(A,M_mitgcm1);

Eigs(1,1) = max(eigen1);
Eigs(1,2) = min(eigen1);%最大最小特征值
Eigs(2,1) = max(eigen2);
Eigs(2,2) = min(eigen2);%最大最小特征值
Eigs(3,1) = max(eigen3);
Eigs(3,2) = min(eigen3);%最大最小特征值
Eigs(4,1) = max(eigen4);
Eigs(4,2) = min(eigen4);%最大最小特征值

[~, iter_mitgcm1] = Pcsi_iter(A,b,x0,M_ones,     tol,Eigs(1,1),Eigs(1,2),iterNum);
[~, iter_mitgcm2] = Pcsi_iter(A,b,x0,M_diagonal, tol,Eigs(2,1),Eigs(2,2),iterNum);
[x1,iter_mitgcm3] = Pcsi_iter(A,b,x0,M_mitgcm,   tol,Eigs(3,1),Eigs(3,2),iterNum);
[x2,iter_mitgcm4] = Pcsi_iter(A,b,x0,M_mitgcm1,  tol,Eigs(4,1),Eigs(4,2),iterNum);


[~, iter_chev1  ] = Chebyshev(A,b,x0,M_ones,     tol,Eigs(1,1),Eigs(1,2),iterNum);
[~, iter_chev2  ] = Chebyshev(A,b,x0,M_diagonal, tol,Eigs(2,1),Eigs(2,2),iterNum);
[~, iter_chev3  ] = Chebyshev(A,b,x0,M_mitgcm,   tol,Eigs(3,1),Eigs(3,2),iterNum);
[~, iter_chev4  ] = Chebyshev(A,b,x0,M_mitgcm1,  tol,Eigs(4,1),Eigs(4,2),iterNum);

x_pcg_diagonal = pcg(-A, -b, tol, 1000);
%[cg2d_x, err_sq, cg2d_m] = mitgcm_cg2d();
[x_cg, iter_cg_Gear ]  = CG_Gear(A, b, x0, M_mitgcm1, tol, iterNum);
% tf = issymmetric(M_diagonal*A);
% d  = eig(A);
% isposdef = all(d>0);
norm(x2-x,2)

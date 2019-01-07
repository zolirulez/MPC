function [Phi,Gamma1,Gamma2] =  delayjob(A,B,m,Ts)
nx = length(A);
nu = size(B,2);
Matrix1 = expm([A B; zeros(nu,nx+nu)]*(1-m)*Ts);
Matrix2 = expm([A B; zeros(nu,nx+nu)]*m*Ts);
Phi = Matrix1(1:nx,1:nx)*Matrix2(1:nx,1:nx);
Gamma1 = Matrix2(1:nx,1:nx)*Matrix1(1:nx,nx+1:end);
Gamma2 = Matrix2(1:nx,nx+1:end);
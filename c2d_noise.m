function [Abar,Qbar]=c2d_noise(A,G,Ts)
[nx,~]=size(G);
M = [-A' G*G'; zeros(nx,nx) A];
Phi = expm(M*Ts);
Abar = Phi(nx+1:nx+nx,nx+1:nx+nx);
Qbar = Abar'*Phi(1:nx,nx+1:nx+nx);
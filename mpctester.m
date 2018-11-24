clearvars
% Kalman Filter 
system.A = [0.9 -1; 0 0.8];
system.B = [0.1 0.1; 2 2];
system.C = eye(2);
system.Cz = [1 1; 1 0];
system.G = [0.1; 2];
noise.R = 0.1*eye(2);
noise.Q = 0.01*eye(1);
noise.S = 0.1*ones(1,2);
initial.x = zeros(2,1);
initial.P = eye(2);
kfType = 'timeinvariant';
kf = KalmanFilter;
horizon = 3;
kf.initialize(system,noise,initial,kfType,horizon)

% ... and then the MPC
horizon = 3;
n.nu = 2;
n.nz = 2;
n.nx = 3;
Gamma = kf.markov;
Uc = zeros(6,1);
% Weights
W.Wz = eye(2);
W.Wu = eye(2);
W.WDu = eye(2);
W.Ws1 = eye(2);
W.Ws2 = eye(2);
W.Wt1 = eye(2);
W.Wt2 = eye(2);
mpc = ModelPredictiveController;
mpc.initialize(W,Gamma,Uc,horizon,n)

DUmin = [-1; -1];
DUmax = [1; 1];
u_1 = zeros(2,1);
b = zeros(horizon*2,1);
R = ones(horizon*2,1);
inputBounds.Umin = -100*ones(horizon*2,1);
inputBounds.Umax = 100*ones(horizon*2,1);
inputBounds.DUmin = -1000*ones(horizon*2,1);
inputBounds.DUmax = 1000*ones(horizon*2,1);
refBounds.Rmin = -100*ones(horizon*2,1);
refBounds.Rmax = 100*ones(horizon*2,1);
optioptions.Algorithm = 'interior-point-convex';
U = mpc.controlCompute(b,R,refBounds,inputBounds,u_1,optioptions)
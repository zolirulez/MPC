%% Exercise 1
% System parameters
K = 1;
tau = 5;
tau1 = tau; tau2 = tau;
beta = 2;
theta = 5;
Kw = 0.1;
tauw = 10;
sigmav = 0.1;

G = K*tf([beta 1],[tau1*tau2 tau1+tau2 1]); % TODO for the delay
H = Kw*tf(1,[tauw 1]);

%% Exercise 2

Ts = 2;
sscu = ss(G);
sscw = ss(H);
Ac = [sscu.A zeros(2,1); zeros(1,2) sscw.A];
Bc = [sscu.B; 0];
Gc = [zeros(2,1); sscw.B];
Czc = [sscu.C sscw.C];
[A,B] = c2d(Ac,Bc,Ts);
[A,G] = c2d(Ac,Gc,Ts);
Cz = Czc;
C = Cz;
% Covariance matrices TODO

%% Exercise 3

system.A = A;
system.B = B;
system.C = C;
system.Cz = Cz;
system.G = G;
% TODO
% noise.R = 0.1*eye(2);
% noise.Q = 0.01*eye(1);
% noise.S = 0.1*ones(1,2);
initial.x = zeros(3,1);
initial.P = eye(3);
kfType = 'stationary';
kf = KalmanFilter;
horizon = 10;
kf.initialize(system,noise,initial,kfType,horizon);
% Print results
kf.Kx
kf.Kw
% Small simulation for showing how to compute x
j = horizon;
u = ones(1,j);
Lr = chol(noise.R);% TODO: the noise variance
y = initial.x + Lr*randn(2,1); 
Q = cell(5,1);
for it = 1:j
    Q{it} = noise.Q;
end
[xf, x1, xj, z] = kf.outputPredictor(u,y,Q,j)

%% Exercise 4
% TODO: THE MARKOV PREDICTOR STILL DOES NOT WORK WELL
kf.markov
kf.obs
kf.obsn
[xf, xj, z] = kf.markovPredictor(u,y)

%% Exercise 5: input constained MPC initialization
% TODO: put the design here
mpc = ModelPredictiveController;
mpc.initialize(W,Gamma,Uc,horizon,n)
%% Exercise 6: compute time steps
% TODO: put the design here
optioptions.Algorithm = 'interior-point-convex';
U = mpc.controlCompute(b,R,refBounds,inputBounds,u_1,optioptions)
%% Exercise 7: closed loop simulation of MPC
% TODO
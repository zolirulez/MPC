clearvars
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

Gtf = K*tf([beta 1],[tau1*tau2 tau1+tau2 1]); % TODO for the delay
Htf = Kw*tf(1,[tauw 1]);

%% Exercise 2

Ts = 2;
sscu = ss(Gtf);
sscw = ss(Htf);
Ac = [sscu.A zeros(2,1); zeros(1,2) sscw.A];
Bc = [sscu.B; 0];
Gc = [zeros(2,1); sscw.B];
Cc = [sscu.C sscw.C];
Czc = [sscu.C 0];
[A,B] = c2d(Ac,Bc,Ts);
[A,G] = c2d(Ac,Gc,Ts);
Cz = Czc;
C = Cc;
% Covariance matrices TODO
Rww = c2d_noise(A,G,Ts);
Rvv = sigmav^2;
Rwv = 0;

%% Exercise 3: Kalman predictor

system.A = A;
system.B = B;
system.C = C;
system.Cz = Cz;
system.G = G;
noise.R = Rvv;
noise.Q = 1;
noise.S = Rwv;
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
Lr = chol(noise.R);
y = C*initial.x + Lr'*randn(1,1); 
Q = cell(5,1);
for it = 1:j
    Q{it} = noise.Q;
end
[xf, x1, xj, z] = kf.outputPredictor(u,y,Q,j)

%% Exercise 4: Markov predictor
kf.Markov
kf.obs
kf.obsn
[xf, xj, zj] = kf.markovPredictor(u,y)

%% Exercise 5: input constained MPC initialization
mpc = ModelPredictiveController;
% Preinitialization
n.nx = 3;
n.nu = 1;
n.nz = 1;
n.no = 2;
n.nscl = 0;
n.nscu = 0;
mpc.preinitialize(horizon,n);
% Initialization
W = cell(2,1);
Ucoeff = cell(2,1);
ExtraFeatures.useInputConstraints = true;
ExtraFeatures.useInputRateConstraints = true;
ExtraFeatures.useSoftOutputConstraints = false;
W{1} = 10;
W{2} = 1;
Ucoeff{1} = kf.Markov;
Ucoeff{2} = mpc.Lambda;
u_1 = 0;
mpc.initialize(W,Ucoeff,ExtraFeatures,u_1)
%% Exercise 6: compute time steps
c = cell(2,1);
b = kf.obs*kf.xf;
Z = kron(ones(10,1),mean(reshape(kf.zj,1,10),2));
c{1} = Z-b;
c{2} = mpc.I0*mpc.u_1;
ExtraFeatures.Bounds.Umin = -100*ones(10,1);
ExtraFeatures.Bounds.Umax = 100*ones(10,1);
ExtraFeatures.Bounds.DUmin = -100*ones(10,1); 
ExtraFeatures.Bounds.DUmax = 100*ones(10,1);
ExtraFeatures.b = b;
optioptions.Algorithm = 'interior-point-convex';
[u, U] = mpc.controlCompute(c,ExtraFeatures,optioptions)
%% Exercise 7: closed loop simulation of MPC
x = zeros(3,1);
Lq = chol(noise.Q); 
itmax = Ts*100;
record.t = 1:Ts:itmax;
record.x = zeros(kf.nx,length(record.t));
record.z = zeros(kf.nz,length(record.t));
for it = 1:length(record.t)
    % Kalman Filter, Markov predictor
    kf.markovPredictor(U,y);
    % Model Predictive Control
    b = kf.obs*kf.xf;
    Z = kron(ones(10,1),mean(reshape(kf.zj,1,10),2));
    c{1} = Z-b;
    c{2} = mpc.I0*mpc.u_1;
    ExtraFeatures.b = kf.obs;
    [u U] = mpc.controlCompute(c,ExtraFeatures,optioptions);
    % Physical system and measurement
    v = Lr'*randn(1,1);
    w = Lq'*randn(1,1);
    x = A*x + B*u + G*w;
    z = Cz*x;
    y = z + v;
    % Saving results
    record.x(1:kf.nx,it) = x;
    record.z(1:kf.nz,it) = z;
end
figure(1)
plot(record.t,record.x)
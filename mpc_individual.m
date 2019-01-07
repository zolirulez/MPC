clearvars
close all
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
Ts = 2;

Gtf = K*tf([beta 1],[tau1*tau2 tau1+tau2 1]); 
Htf = Kw*tf(1,[tauw 1]);

%% Exercise 2
m = theta/Ts - floor(theta/Ts);
sscu = ss(Gtf);
sscw = ss(Htf);
% Augmentation with noise
Ac = [sscu.A zeros(2,1); zeros(1,2) sscw.A];
Bc = [sscu.B; 0];
Gc = [zeros(2,1); sscw.B];
Czc = [sscu.C sscw.C];
Dc = sscu.D;
Cz = [Czc Dc zeros(1,2)];
C = [Czc Dc zeros(1,2)];
% Delay calculations
[Phi,Gamma1,Gamma2] =  delayjob(Ac,Bc,m,Ts);
A = [Phi Gamma1 Gamma2 zeros(3,1); zeros(3,4) eye(3,2)];
B = [zeros(5,1); 1];
% Covariance matrices
[~,Rww] = c2d_noise(Ac,Gc,Ts);
Rww = [Rww zeros(3,3); zeros(3,6)];
Rvv = sigmav^2;
Rwv = zeros(6,1);
G = eye(6);
% Reporting matrices
syms dummy
disp('Exercise2')
latex(vpa(subs(A+dummy,dummy,0),4))
latex(vpa(subs(B+dummy,dummy,0),4))
latex(vpa(subs(G+dummy,dummy,0),4))
latex(vpa(subs(Cz+dummy,dummy,0),4))
latex(vpa(subs(Rww+dummy,dummy,0),4))
latex(vpa(subs(Rwv+dummy,dummy,0),4))
latex(vpa(subs(Rvv+dummy,dummy,0),4))

%% Exercise 3: Kalman predictor

system.A = A;
system.B = B;
system.C = C;
system.Cz = Cz;
system.G = G;
noise.R = Rvv;
noise.Q = Rww;
noise.S = Rwv;
initial.x = zeros(6,1);
initial.P = eye(6);
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
Lr = chol(noise.R,'lower');
y = C*initial.x + Lr'*randn(1,1); 
Q = cell(j,1);
for it = 1:j
    Q{it} = noise.Q;
end
[xf, x1, xj, z] = kf.outputPredictor(u,y,Q,j);
% Reporting matrices
disp('Exercise3')
latex(vpa(subs(kf.Kx+dummy,dummy,0),4))
latex(vpa(subs(kf.Kw+dummy,dummy,0),4))

%% Exercise 4: Markov predictor
% Report matrices
disp('Exercise4')
latex(vpa(subs(kf.obs+dummy,dummy,0),4))
latex(vpa(subs(kf.obsn+dummy,dummy,0),4))
latex(vpa(subs(kf.Markov+dummy,dummy,0),4))
[xf, xj, zj] = kf.markovPredictor(u,y);

%% Exercise 5: input constained MPC initialization
mpc = ModelPredictiveController;
% Preinitialization
n = kf.stateSpaceDimensions;
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
W{2} = 5;
Ucoeff{1} = kf.Markov;
Ucoeff{2} = mpc.Lambda;
u_1 = 0;
optioptions.Algorithm = 'interior-point-convex';
optioptions.Display = 'off';
mpc.initialize(W,Ucoeff,n,ExtraFeatures,u_1,optioptions);
%% Exercise 6: compute time steps
c = cell(2,1);
Z = ones(10,1);
c{1} = Z - kf.b;
c{2} = mpc.I0*mpc.u_1;
ExtraFeatures.Bounds.Umin = -2*ones(10,1);
ExtraFeatures.Bounds.Umax = 2*ones(10,1);
ExtraFeatures.Bounds.DUmin = -10*ones(10,1); 
ExtraFeatures.Bounds.DUmax = 10*ones(10,1);
ExtraFeatures.b = kf.b;
[u, U] = mpc.controlCompute(c,ExtraFeatures);
% Report matrices
disp('Exercise5')
latex(vpa(subs(mpc.H+dummy,dummy,0),4))
latex(vpa(subs(-mpc.M{1}+dummy,dummy,0),4))
latex(vpa(subs(mpc.M{1}+dummy,dummy,0),4))
latex(vpa(subs(mpc.M{2}+dummy,dummy,0),4))
latex(vpa(subs(mpc.l+dummy,dummy,0),4))
latex(vpa(subs(mpc.u+dummy,dummy,0),4))
latex(vpa(subs(mpc.A+dummy,dummy,0),4))
latex(vpa(subs(mpc.I0+dummy,dummy,0),4))
latex(vpa(subs(ExtraFeatures.Bounds.DUmin+dummy,dummy,0),4))
latex(vpa(subs(ExtraFeatures.Bounds.DUmax+dummy,dummy,0),4))
disp('Exercise6')
%% Exercise 7: closed loop simulation of MPC
x = zeros(6,1);
Lq = diag([0 0 sqrt(noise.Q(3,3)) 0 0 0]);
itmax = Ts*100;
R = ones(1,itmax/Ts+horizon);
for it = 1:(itmax/Ts+horizon)
    if rem(floor(it/(itmax/10)),2)
        R(:,it) = -R(:,it);
    end
end
record.t = 0:Ts:itmax-Ts+horizon*Ts;
record.x = zeros(kf.nx,length(record.t));
record.z = zeros(kf.nz,length(record.t));
record.u = zeros(kf.nu,length(record.t));
fighand = figure;
clf
subplot(211)
plot(record.t,R,'b--')
axis([0 itmax+horizon*Ts+Ts -2 2])
hand11 = animatedline('Color','k');
hand12 = animatedline('Color','r');
subplot(212)
plot(record.t,diag([2 -2])*ones(2,itmax/Ts+horizon),'r--')
axis([0 itmax+horizon*Ts+Ts -3 3])
xlabel('Time [-]')
hand2 = animatedline('Color','k');
% Initial values
u = 0;
z = 0;
y = 0;
kf.reinitialize(initial,horizon)
U = zeros(j,1);
% Seed
Seed = 12345;
randn('state',Seed);
for it = 1:itmax/Ts
    % Kalman Filter, Markov predictor
    kf.markovPredictor(U,y);
    % Model Predictive Control
    Z = R(1,it:it+horizon-1)';
    c{1} = Z-kf.b;
    c{2} = mpc.I0*mpc.u_1;
    ExtraFeatures.b = kf.b;
    [u, U] = mpc.controlCompute(c,ExtraFeatures);
    % Print prediction
    subplot(211)
    hold on
    plot((it+1)*Ts-2*Ts:Ts:(it+kf.j)*Ts-2*Ts,kf.zj,'g-')
    hold off
    subplot(212)
    hold on
    plot((it+1)*Ts-2*Ts:Ts:(it+kf.j)*Ts-2*Ts,U,'g-')
    hold off
    % Animation
    pause(0.01)
    addpoints(hand11,record.t(it),z)
    addpoints(hand12,record.t(it),y)
    addpoints(hand2,record.t(it),u)
    % Physical system and measurement
    v = Lr*randn(1,1);
    w = Lq*randn(6,1);
    x = A*x + B*u + G*w;
    z = C*x;
    y = z + v;
    % Saving results
    record.x(1:kf.nx,it) = x;
    record.z(1:kf.nz,it) = z;
    record.u(1:kf.nz,it) = u;
end
addpoints(hand11,record.t(it),z)
addpoints(hand12,record.t(it),y)
addpoints(hand2,record.t(it),u)
subplot(211)
legend('Reference r','Output z','Measurement y','Markov prediction Z')
subplot(212)
legend('Lower bound','Upper bound','Input u','Future input U')
uistack(hand12,'top');
uistack(hand11,'top');
uistack(hand2,'top');
saveas(fighand,[pwd '\figures\individual'],'png')